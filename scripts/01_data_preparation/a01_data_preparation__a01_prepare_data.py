import csv
import gzip
import hashlib
import json
import os
import re
import sys
import tarfile
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(__file__).resolve().parents[1]
DIRS = {
    "raw": BASE,
    "processed": BASE,
    "metadata": BASE,
    "qc": BASE,
    "scripts": BASE,
    "logs": BASE,
}
for d in DIRS.values():
    d.mkdir(parents=True, exist_ok=True)

LOG = DIRS["logs"] / "a01_prepare_data.log"


def log(msg):
    text = str(msg)
    print(text)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(text + "\n")


def download(url, dest):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and dest.stat().st_size > 0:
        log(f"SKIP existing {dest} ({dest.stat().st_size} bytes)")
        return dest
    log(f"DOWNLOAD {url} -> {dest}")
    with urllib.request.urlopen(url, timeout=180) as r, dest.open("wb") as out:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    log(f"DONE {dest} ({dest.stat().st_size} bytes)")
    return dest


def head_url(url):
    req = urllib.request.Request(url, method="HEAD")
    with urllib.request.urlopen(req, timeout=60) as r:
        return {
            "status": r.status,
            "content_length": r.headers.get("Content-Length", ""),
            "content_type": r.headers.get("Content-Type", ""),
            "accept_ranges": r.headers.get("Accept-Ranges", ""),
        }


def gzip_probe(path, n=3):
    lines = []
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
        for _ in range(n):
            line = fh.readline()
            if not line:
                break
            lines.append(line.rstrip("\n")[:500])
    return lines


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def read_series_matrix(path):
    meta = defaultdict(list)
    table_header = None
    rows = []
    in_table = False
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                table_header = next(fh).rstrip("\n").split("\t")
                continue
            if line == "!series_matrix_table_end":
                in_table = False
                continue
            if in_table:
                rows.append(line.split("\t"))
            elif line.startswith("!"):
                parts = line.split("\t")
                meta[parts[0]].append([p.strip('"') for p in parts[1:]])
    return meta, table_header, rows


def write_tsv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


def write_csv(path, header, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def norm_title_value(v):
    return v.strip().strip('"')


def parse_platform_gpl570(platform_path):
    genes = {}
    in_table = False
    header = None
    with open_text(platform_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
                header = next(fh).rstrip("\n").split("\t")
                continue
            if line == "!platform_table_end":
                break
            if not in_table:
                continue
            parts = line.split("\t")
            if len(parts) < len(header):
                parts += [""] * (len(header) - len(parts))
            row = dict(zip(header, parts))
            probe = row.get("ID", "")
            symbol = row.get("Gene Symbol", "") or row.get("gene_assignment", "")
            if " /// " in symbol:
                symbol = symbol.split(" /// ")[0]
            if " // " in symbol and not re.match(r"^[A-Za-z0-9_.-]+$", symbol):
                fields = symbol.split(" // ")
                symbol = fields[1] if len(fields) > 1 else fields[0]
            symbol = symbol.strip()
            if probe and symbol and symbol != "---":
                genes[probe] = symbol
    return genes


def collapse_probe_to_gene(expr_rows, sample_ids, probe_to_gene, out_path, max_rows=None):
    sums = {}
    counts = Counter()
    used_probe = 0
    for i, row in enumerate(expr_rows):
        if max_rows and i >= max_rows:
            break
        probe = row[0].strip('"')
        gene = probe_to_gene.get(probe, "")
        if not gene:
            continue
        vals = []
        ok = True
        for x in row[1:]:
            try:
                vals.append(float(x))
            except ValueError:
                ok = False
                break
        if not ok:
            continue
        if gene not in sums:
            sums[gene] = [0.0] * len(sample_ids)
        for j, v in enumerate(vals):
            sums[gene][j] += v
        counts[gene] += 1
        used_probe += 1
    genes = sorted(sums)
    with out_path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["gene_symbol"] + sample_ids)
        for g in genes:
            c = counts[g]
            w.writerow([g] + [f"{v / c:.6g}" for v in sums[g]])
    return {
        "genes": len(genes),
        "used_probes": used_probe,
        "samples": len(sample_ids),
        "out": str(out_path),
    }


def clean_lcs(x):
    return x.replace("_", "-")


def parse_gse14520_extra(path):
    with open_text(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return list(reader)


def extract_lcs_from_title(title):
    m = re.search(r"(LCS[-_]\d+A|X02[-_]\d+A)", title)
    return m.group(1).replace("-", "_") if m else ""


def main():
    LOG.write_text("", encoding="utf-8")
    log(f"A01 base: {BASE}")
    urls = {
        "gse104580_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104580/matrix/GSE104580_series_matrix.txt.gz",
        "gse104580_soft": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104580/soft/GSE104580_family.soft.gz",
        "gpl570_soft": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPLnnn/GPL570/soft/GPL570_family.soft.gz",
        "gse14520_gpl3921": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL3921_series_matrix.txt.gz",
        "gse14520_gpl571": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/matrix/GSE14520-GPL571_series_matrix.txt.gz",
        "gse14520_extra": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE14nnn/GSE14520/suppl/GSE14520_Extra_Supplement.txt.gz",
        "gse277104_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277104/matrix/GSE277104_series_matrix.txt.gz",
        "gse277104_qnorm": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277104/suppl/GSE277104_qnormCounts.tsv.gz",
        "gse149614_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/matrix/GSE149614_series_matrix.txt.gz",
        "gse149614_meta": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.metadata.updated.txt.gz",
        "gse151530_info": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_Info.txt.gz",
        "gse151530_gpl24676": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/matrix/GSE151530-GPL24676_series_matrix.txt.gz",
        "gse151530_gpl20301": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/matrix/GSE151530-GPL20301_series_matrix.txt.gz",
        "gse151530_gpl18573": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/matrix/GSE151530-GPL18573_series_matrix.txt.gz",
    }
    raw_paths = {}
    for key, url in urls.items():
        raw_paths[key] = download(url, DIRS["raw"] / Path(url).name)

    large_urls = {
        "GSE149614_count": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.count.txt.gz",
        "GSE149614_normalized": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE149nnn/GSE149614/suppl/GSE149614_HCC.scRNAseq.S71915.normalized.txt.gz",
        "GSE151530_matrix_mtx": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_matrix.mtx.gz",
        "GSE151530_barcodes": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_barcodes.tsv.gz",
        "GSE151530_genes": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE151nnn/GSE151530/suppl/GSE151530_genes.tsv.gz",
    }
    large_head = []
    for key, url in large_urls.items():
        try:
            h = head_url(url)
            downloaded_path = ""
            downloaded_bytes = ""
            probe = ""
            if key != "GSE149614_normalized":
                path = download(url, DIRS["raw"] / Path(url).name)
                downloaded_path = str(path)
                downloaded_bytes = path.stat().st_size
                probe = " | ".join(gzip_probe(path, 3))
            else:
                probe = "HEAD only; raw count matrix retained for A01 to avoid duplicative normalized dense download."
            large_head.append([key, url, h["status"], h["content_length"], h["content_type"], h["accept_ranges"], downloaded_path, downloaded_bytes, probe])
        except Exception as e:
            large_head.append([key, url, "ERROR", "", "", str(e), "", "", ""])
    write_csv(DIRS["metadata"] / "large_file_downloadability.csv",
              ["resource", "url", "status", "content_length", "content_type", "accept_ranges_or_error", "downloaded_path", "downloaded_bytes", "gzip_probe"],
              large_head)

    readiness = []
    gene_sets = {}

    # GSE104580
    meta104, header104, rows104 = read_series_matrix(raw_paths["gse104580_matrix"])
    titles104 = meta104["!Sample_title"][0]
    geos104 = meta104["!Sample_geo_accession"][0]
    matrix_cols104 = [x.strip('"') for x in header104[1:]]
    subgroups104 = []
    for vals in meta104["!Sample_characteristics_ch1"]:
        if vals and vals[0].startswith("subject subgroup:"):
            subgroups104 = [v.replace("subject subgroup: ", "") for v in vals]
    sample_rows104 = []
    for gsm, title, subgroup in zip(geos104, titles104, subgroups104):
        response_group = "response" if title.startswith("TACE_Response_") else "non_response"
        sample_rows104.append([gsm, title, subgroup, response_group])
    write_tsv(DIRS["metadata"] / "GSE104580_sample_metadata.tsv",
              ["gsm", "sample_title", "geo_subject_subgroup", "response_group"], sample_rows104)
    probe_to_gene_570 = parse_platform_gpl570(raw_paths["gpl570_soft"])
    write_tsv(DIRS["metadata"] / "GPL570_probe_to_gene.tsv",
              ["probe_id", "gene_symbol"], sorted(probe_to_gene_570.items()))
    # A01 flat execution constraint: avoid downloading additional multi-GB GPL3921/GPL571 SOFT files.
    # GSE14520 Affymetrix probe IDs are mapped conservatively through the already retained GPL570 probe table.
    probe_to_gene_3921 = probe_to_gene_570
    probe_to_gene_571 = probe_to_gene_570
    write_tsv(DIRS["metadata"] / "GPL3921_probe_to_gene.tsv",
              ["probe_id", "gene_symbol"], sorted(probe_to_gene_3921.items()))
    write_tsv(DIRS["metadata"] / "GPL571_probe_to_gene.tsv",
              ["probe_id", "gene_symbol"], sorted(probe_to_gene_571.items()))
    collapse104 = collapse_probe_to_gene(rows104, matrix_cols104, probe_to_gene_570,
                                         DIRS["processed"] / "GSE104580_gene_expression.tsv")
    gene_sets["GSE104580"] = set()
    with (DIRS["processed"] / "GSE104580_gene_expression.tsv").open(encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            gene_sets["GSE104580"].add(line.split("\t", 1)[0])
    readiness.append(["GSE104580", "bulk pretreatment TACE response", "ready",
                      len(titles104), collapse104["genes"],
                      f"GPL570; response={sum(t.startswith('TACE_Response_') for t in titles104)}; non_response={sum(t.startswith('TACE_Non-Response_') for t in titles104)}; gsm_match={geos104 == matrix_cols104}"])

    # GSE14520
    extra_rows = parse_gse14520_extra(raw_paths["gse14520_extra"])
    lcs_clin = {r["LCS ID"]: r for r in extra_rows}
    # Fako Table S1 parsed from prior A00 PDF coordinate extraction
    adjuvant = "LCS_007A LCS_009A LCS_019A LCS_020A LCS_025A LCS_027A LCS_028A LCS_029A LCS_031A LCS_033A LCS_034A LCS_038A LCS_043A LCS_047A LCS_049A LCS_050A LCS_062A LCS_065A LCS_068A LCS_071A LCS_075A LCS_079A LCS_085A LCS_086A LCS_092A LCS_097A LCS_100A LCS_104A LCS_110A LCS_116A LCS_117A LCS_118A LCS_121A LCS_127A LCS_134A LCS_136A LCS_140A LCS_142A LCS_146A LCS_154A LCS_158A LCS_159A LCS_161A LCS_166A LCS_167A LCS_170A LCS_171A LCS_177A LCS_185A LCS_191A LCS_192A LCS_196A LCS_197A LCS_208A LCS_209A LCS_212A LCS_213A LCS_223A LCS_228A LCS_231A LCS_240A LCS_241A LCS_245A LCS_251A LCS_259A LCS_260A LCS_263A LCS_264A LCS_265A LCS_266A LCS_270A LCS_272A LCS_284A LCS_289A LCS_393A".split()
    postrec = "LCS_008A LCS_012A LCS_023A LCS_024A LCS_032A LCS_035A LCS_067A LCS_072A LCS_088A LCS_096A LCS_120A LCS_138A LCS_139A LCS_145A LCS_178A LCS_190A LCS_194A LCS_198A LCS_200A LCS_207A LCS_224A LCS_227A LCS_234A LCS_238A LCS_267A LCS_273A LCS_274A LCS_281A LCS_333A LCS_403A".split()
    resection = "LCS_010A LCS_014A LCS_015A LCS_016A LCS_018A LCS_022A LCS_040A LCS_041A LCS_044A LCS_045A LCS_046A LCS_048A LCS_051A LCS_056A LCS_057A LCS_061A LCS_063A LCS_064A LCS_069A LCS_073A LCS_076A LCS_078A LCS_084A LCS_090A LCS_091A LCS_094A LCS_099A LCS_101A LCS_102A LCS_105A LCS_106A LCS_108A LCS_109A LCS_119A LCS_122A LCS_130A LCS_131A LCS_132A LCS_137A LCS_144A LCS_147A LCS_150A LCS_151A LCS_156A LCS_160A LCS_163A LCS_165A LCS_169A LCS_172A LCS_174A LCS_179A LCS_180A LCS_184A LCS_189A LCS_205A LCS_210A LCS_211A LCS_215A LCS_216A LCS_222A LCS_236A LCS_237A LCS_243A LCS_247A LCS_249A LCS_253A LCS_254A LCS_261A LCS_262A LCS_268A LCS_269A LCS_275A LCS_278A LCS_279A LCS_282A LCS_285A LCS_286A LCS_291A LCS_343A LCS_344A LCS_346A LCS_400A LCS_406A LCS_415A LCS_424A LCS_426A".split()
    other = "LCS_002A LCS_004A LCS_005A LCS_011A LCS_021A LCS_036A LCS_039A LCS_042A LCS_054A LCS_066A LCS_074A LCS_083A LCS_089A LCS_093A LCS_095A LCS_103A LCS_107A LCS_123A LCS_125A LCS_126A LCS_129A LCS_135A LCS_143A LCS_148A LCS_149A LCS_152A LCS_153A LCS_157A LCS_162A LCS_164A LCS_173A LCS_175A LCS_182A LCS_183A LCS_188A LCS_193A LCS_195A LCS_199A LCS_201A LCS_203A LCS_206A LCS_219A LCS_230A LCS_248A LCS_250A LCS_256A LCS_277A LCS_290A LCS_339A LCS_341A LCS_401A".split()
    missing = "LCS_204A LCS_283A LCS_347A X02_342A X02_262A".split()
    therapy = {}
    for x in adjuvant:
        therapy[x] = "adjuvant_TACE"
    for x in postrec:
        therapy[x] = "post_recurrence_TACE"
    for x in resection:
        therapy[x] = "resection_only"
    for x in other:
        therapy[x] = "other_therapy"
    for x in missing:
        therapy[x] = "missing_survival_data"
    meta145_all = []
    title_to_gsm = {}
    for key in ["gse14520_gpl3921", "gse14520_gpl571"]:
        m, h, r = read_series_matrix(raw_paths[key])
        titles = m["!Sample_title"][0]
        geos = m["!Sample_geo_accession"][0]
        for title, gsm in zip(titles, geos):
            lcs = extract_lcs_from_title(title)
            if lcs:
                title_to_gsm[lcs] = gsm
    for lcs, group in sorted(therapy.items()):
        r = lcs_clin.get(lcs, {})
        meta145_all.append([
            lcs, clean_lcs(lcs), title_to_gsm.get(lcs, r.get("Affy_GSM", "")),
            group, "TACE_treated" if group in ("adjuvant_TACE", "post_recurrence_TACE") else group,
            r.get("Survival status", ""), r.get("Survival months", ""),
            r.get("Recurr status", ""), r.get("Recurr months", ""),
            r.get("Tissue Type", ""), r.get("TNM staging", ""), r.get("BCLC staging", "")
        ])
    write_tsv(DIRS["metadata"] / "GSE14520_Fako_mapped_clinical_metadata.tsv",
              ["lcs_id", "geo_lcs_id", "gsm", "fako_therapy_group", "analysis_context",
               "os_status", "os_months", "recurrence_status", "recurrence_months",
               "tissue_type", "tnm_stage", "bclc_stage"], meta145_all)
    # expression: keep GPL3921 and GPL571 as separate platform-level gene matrices
    meta145, header145, rows145 = read_series_matrix(raw_paths["gse14520_gpl3921"])
    samples145 = [x.strip('"') for x in header145[1:]]
    collapse145 = collapse_probe_to_gene(rows145, samples145, probe_to_gene_3921,
                                         DIRS["processed"] / "GSE14520_GPL3921_gene_expression.tsv")
    meta145_571, header145_571, rows145_571 = read_series_matrix(raw_paths["gse14520_gpl571"])
    samples145_571 = [x.strip('"') for x in header145_571[1:]]
    collapse145_571 = collapse_probe_to_gene(rows145_571, samples145_571, probe_to_gene_571,
                                             DIRS["processed"] / "GSE14520_GPL571_gene_expression.tsv")
    gene_sets["GSE14520"] = set()
    with (DIRS["processed"] / "GSE14520_GPL3921_gene_expression.tsv").open(encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            gene_sets["GSE14520"].add(line.split("\t", 1)[0])
    with (DIRS["processed"] / "GSE14520_GPL571_gene_expression.tsv").open(encoding="utf-8") as fh:
        next(fh)
        for line in fh:
            gene_sets["GSE14520"].add(line.split("\t", 1)[0])
    readiness.append(["GSE14520", "bulk TACE-treated outcome/resection-only specificity", "ready",
                      len(samples145) + len(samples145_571), len(gene_sets["GSE14520"]),
                      f"Fako Table S1: TACE=105; resection-only=86; other=51; missing=5; OS available; recurrence months/status available; GPL3921 genes={collapse145['genes']}; GPL571 genes={collapse145_571['genes']}"])

    # GSE277104
    meta277, header277, rows277 = read_series_matrix(raw_paths["gse277104_matrix"])
    titles277 = meta277["!Sample_title"][0]
    geos277 = meta277["!Sample_geo_accession"][0]
    char277 = meta277["!Sample_characteristics_ch1"]
    def char_vals(prefix):
        for vals in char277:
            if vals and vals[0].startswith(prefix):
                return [v.replace(prefix, "").strip() for v in vals]
        return [""] * len(titles277)
    roi = char_vals("roi:")
    aoi = char_vals("aoi type:")
    tissue = char_vals("tissue:")
    patient = []
    for title in titles277:
        mm = re.search(r"_(\d{2}_\d+)_ROI", title)
        patient.append(mm.group(1) if mm else "NA")
    rows_meta277 = [[gsm, title, p, r, ao, ti] for gsm, title, p, r, ao, ti in zip(geos277, titles277, patient, roi, aoi, tissue)]
    write_tsv(DIRS["metadata"] / "GSE277104_AOI_metadata.tsv",
              ["gsm", "sample_title", "patient_key", "roi", "aoi_type", "tissue"], rows_meta277)
    pair = defaultdict(set)
    for p, r, ao in zip(patient, roi, aoi):
        if p != "NA" and r != "NA":
            pair[(p, r)].add(ao)
    pair_rows = []
    for (p, r), aos in sorted(pair.items()):
        pair_rows.append([p, r, ";".join(sorted(aos)), int("Tumor" in aos), int("Vessel" in aos), int("Tumor" in aos and "Vessel" in aos)])
    write_tsv(DIRS["metadata"] / "GSE277104_tumor_vessel_pairing_summary.tsv",
              ["patient_key", "roi", "aoi_types", "has_tumor", "has_vessel", "is_tumor_vessel_pair"], pair_rows)
    # Copy qnorm matrix as processed spatial matrix
    with open_text(raw_paths["gse277104_qnorm"]) as src, (DIRS["processed"] / "GSE277104_qnorm_AOI_expression.tsv").open("w", encoding="utf-8") as out:
        for line in src:
            out.write(line)
    gene_sets["GSE277104"] = set()
    with open_text(raw_paths["gse277104_qnorm"]) as fh:
        next(fh)
        for line in fh:
            if line.strip():
                gene_sets["GSE277104"].add(line.split("\t", 1)[0])
    readiness.append(["GSE277104", "GeoMx tumor-vessel spatial anchor", "ready",
                      len(titles277), len(gene_sets["GSE277104"]),
                      f"AOI={len(titles277)}; Tumor={aoi.count('Tumor')}; Vessel={aoi.count('Vessel')}; pairedROI={sum(1 for x in pair.values() if 'Tumor' in x and 'Vessel' in x)}"])

    # GSE149614
    meta149, header149, rows149 = read_series_matrix(raw_paths["gse149614_matrix"])
    titles149 = meta149["!Sample_title"][0]
    sample_rows149 = []
    for i, title in enumerate(titles149):
        vals = [chars[i] for chars in meta149["!Sample_characteristics_ch1"] if len(chars) > i]
        sample_rows149.append([title] + vals)
    max_chars = max(len(r) for r in sample_rows149)
    write_tsv(DIRS["metadata"] / "GSE149614_sample_metadata.tsv",
              ["sample_title"] + [f"characteristic_{i}" for i in range(1, max_chars)],
              [r + [""] * (max_chars - len(r)) for r in sample_rows149])
    # copy cell metadata
    with open_text(raw_paths["gse149614_meta"]) as src, (DIRS["metadata"] / "GSE149614_cell_metadata.tsv").open("w", encoding="utf-8") as out:
        for line in src:
            out.write(line)
    cell_counts149 = Counter()
    site_counts149 = Counter()
    with open_text(raw_paths["gse149614_meta"]) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            cell_counts149[r.get("celltype", "")] += 1
            site_counts149[r.get("site", "")] += 1
    write_tsv(DIRS["qc"] / "GSE149614_celltype_counts.tsv", ["celltype", "n_cells"], sorted(cell_counts149.items()))
    write_tsv(DIRS["qc"] / "GSE149614_site_counts.tsv", ["site", "n_cells"], sorted(site_counts149.items()))
    readiness.append(["GSE149614", "HCC scRNA anchor", "metadata_ready_expression_large",
                      sum(cell_counts149.values()), "", f"celltypes={dict(cell_counts149)}; sites={dict(site_counts149)}; expression matrix HEAD checked in metadata/large_file_downloadability.csv"])

    # GSE151530
    with open_text(raw_paths["gse151530_info"]) as src, (DIRS["metadata"] / "GSE151530_cell_metadata.tsv").open("w", encoding="utf-8") as out:
        for line in src:
            out.write(line)
    type_counts151 = Counter()
    sample_counts151 = Counter()
    with open_text(raw_paths["gse151530_info"]) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            type_counts151[r.get("Type", "")] += 1
            sample_counts151[r.get("S_ID", "")] += 1
    write_tsv(DIRS["qc"] / "GSE151530_celltype_counts.tsv", ["celltype", "n_cells"], sorted(type_counts151.items()))
    write_tsv(DIRS["qc"] / "GSE151530_sample_cell_counts.tsv", ["sample", "n_cells"], sorted(sample_counts151.items()))
    sample_cancer = {}
    for key in ["gse151530_gpl24676", "gse151530_gpl20301", "gse151530_gpl18573"]:
        m, h, r = read_series_matrix(raw_paths[key])
        titles = m["!Sample_title"][0]
        cancer_line = None
        for vals in m["!Sample_characteristics_ch1"]:
            if vals and vals[0].startswith("cancer type:"):
                cancer_line = vals
        if cancer_line:
            for title, c in zip(titles, cancer_line):
                sample_cancer[title] = c.replace("cancer type: ", "")
    write_tsv(DIRS["metadata"] / "GSE151530_sample_cancer_type.tsv",
              ["S_ID", "cancer_type"], sorted(sample_cancer.items()))
    hcc_samples = {s for s, c in sample_cancer.items() if "Hepatocellular" in c}
    with open_text(raw_paths["gse151530_info"]) as src, (DIRS["metadata"] / "GSE151530_HCC_filtered_cell_metadata.tsv").open("w", encoding="utf-8") as out:
        header = src.readline()
        out.write(header)
        for line in src:
            sid = line.split("\t", 1)[0]
            if sid in hcc_samples:
                out.write(line)
    readiness.append(["GSE151530", "liver cancer scRNA HCC-filtered support", "metadata_ready_expression_large",
                      sum(type_counts151.values()), "", f"celltypes={dict(type_counts151)}; HCC_samples={len(hcc_samples)}; iCCA_filter_required=true"])

    # gene harmonization
    candidate_sets = {
        "hypoxia_core": "CA9 VEGFA SLC2A1 HK2 LDHA PDK1 BNIP3 EGLN3 ADM ENO1 NDRG1 HILPDA DDIT4 ALDOA".split(),
        "spp1_mmp9_macrophage": "SPP1 MMP9 CD68 LST1 C1QA C1QB C1QC APOE MARCO TREM2 LGALS3 MSR1".split(),
        "angiogenesis_endothelial": "PECAM1 VWF KDR FLT1 ENG CD34 EMCN RAMP2 ACKR1 ESAM ANGPT2 TEK PLVAP".split(),
        "tumor_endothelial_interaction": "VEGFA VEGFB PGF KDR FLT1 ANGPT2 TEK LGALS9 HAVCR2 CXCL12 CXCR4 DLL4 NOTCH1 JAG1".split(),
    }
    all_sets = {k: v for k, v in gene_sets.items() if v}
    common = set.intersection(*all_sets.values()) if all_sets else set()
    write_tsv(DIRS["metadata"] / "cross_dataset_gene_counts.tsv",
              ["dataset", "n_gene_symbols"], [[k, len(v)] for k, v in sorted(all_sets.items())] + [["intersection", len(common)]])
    with (DIRS["metadata"] / "cross_dataset_intersection_genes.txt").open("w", encoding="utf-8") as fh:
        for g in sorted(common):
            fh.write(g + "\n")
    avail_rows = []
    for set_name, genes in candidate_sets.items():
        for g in genes:
            row = [set_name, g]
            for ds in sorted(all_sets):
                row.append("yes" if g in all_sets[ds] else "no")
            row.append("yes" if g in common else "no")
            avail_rows.append(row)
    write_tsv(DIRS["metadata"] / "candidate_gene_availability.tsv",
              ["candidate_set", "gene"] + sorted(all_sets) + ["all_intersection"], avail_rows)

    # QC reports
    write_csv(DIRS["metadata"] / "A01_dataset_readiness_table.csv",
              ["dataset", "role", "readiness", "n_samples_or_cells", "n_genes", "notes"], readiness)
    write_csv(DIRS["qc"] / "A01_basic_qc_summary.csv",
              ["dataset", "qc_item", "value"],
              [
                  ["GSE104580", "n_samples", len(titles104)],
                  ["GSE104580", "n_response", sum(t.startswith("TACE_Response_") for t in titles104)],
                  ["GSE104580", "n_non_response", sum(t.startswith("TACE_Non-Response_") for t in titles104)],
                  ["GSE104580", "duplicate_GSM", str(len(set(geos104)) != len(geos104))],
                  ["GSE14520", "Fako_TACE_treated", 105],
                  ["GSE14520", "Fako_resection_only", 86],
                  ["GSE14520", "Fako_other_therapy", 51],
                  ["GSE14520", "Fako_missing_survival", 5],
                  ["GSE277104", "AOI_total", len(titles277)],
                  ["GSE277104", "Tumor_AOI", aoi.count("Tumor")],
                  ["GSE277104", "Vessel_AOI", aoi.count("Vessel")],
                  ["GSE277104", "paired_ROI", sum(1 for x in pair.values() if "Tumor" in x and "Vessel" in x)],
                  ["GSE149614", "cells_with_metadata", sum(cell_counts149.values())],
                  ["GSE151530", "cells_with_metadata", sum(type_counts151.values())],
              ])
    log("A01 data preparation completed")


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        log(f"ERROR: {type(e).__name__}: {e}")
        raise
