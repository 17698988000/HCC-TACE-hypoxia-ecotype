import csv
import gzip
import math
import re
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
A02 = ROOT / "results" / "a02_ecotype_construction"
A05 = ROOT / "results" / "a05_secondary_externalization"

COMPONENTS = [
    "malignant_hypoxia_adaptive",
    "endothelial_angiogenesis_interaction",
    "spp1_mmp9_tam_matrix_angiogenic_support",
]
PRIMARY_SCORE = "total_ecotype_equal_component_score"


URLS = {
    "GSE76427_series_matrix": "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz",
    "GPL10558_soft": "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/soft/GPL10558_family.soft.gz",
    "TCGA_LIHC_HiSeqV2": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap/HiSeqV2.gz",
    "TCGA_LIHC_clinical": "https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.LIHC.sampleMap/LIHC_clinicalMatrix",
    "TCGA_LIHC_survival": "https://gdc-hub.s3.us-east-1.amazonaws.com/download/TCGA-LIHC.survival.tsv.gz",
}


def download(url, dest):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    with urllib.request.urlopen(url, timeout=180) as r, dest.open("wb") as out:
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
    return dest


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def write_csv(path, header, rows):
    with (A05 / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def split_line(line):
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]


def read_series_matrix(path):
    meta = defaultdict(list)
    header = None
    rows = []
    in_table = False
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                in_table = True
                continue
            if line == "!series_matrix_table_end":
                break
            parts = split_line(line)
            if in_table:
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
            elif line.startswith("!"):
                meta[parts[0]].append(parts[1:])
    return meta, header, rows


def parse_gpl(path):
    in_table = False
    header = None
    probe_to_gene = {}
    with open_text(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "!platform_table_begin":
                in_table = True
                continue
            if line == "!platform_table_end":
                break
            if not in_table:
                continue
            parts = line.split("\t")
            if header is None:
                header = parts
                continue
            id_idx = header.index("ID") if "ID" in header else None
            gene_idx = None
            for cand in ["Gene Symbol", "Symbol", "ILMN_Gene", "Gene", "gene_assignment"]:
                if cand in header:
                    gene_idx = header.index(cand)
                    break
            if id_idx is None or gene_idx is None or len(parts) <= max(id_idx, gene_idx):
                continue
            raw = parts[gene_idx]
            gene = re.split(r"///|//|;|\s+", raw.strip())[0]
            if gene and gene not in {"---", "NA", "nan"}:
                probe_to_gene[parts[id_idx]] = gene
    return probe_to_gene


def collapse_probe_rows(rows, samples, probe_to_gene):
    sums = {}
    counts = Counter()
    for row in rows:
        probe = row[0]
        gene = probe_to_gene.get(probe)
        if not gene:
            continue
        try:
            vals = np.array([float(x) for x in row[1:]], dtype=float)
        except ValueError:
            continue
        if gene not in sums:
            sums[gene] = np.zeros(len(samples), dtype=float)
        sums[gene] += vals
        counts[gene] += 1
    return {g: sums[g] / counts[g] for g in sums}


def load_tcga_hiseq(path):
    expr = {}
    with open_text(path) as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        for row in reader:
            if not row:
                continue
            gene = row[0].split("|")[0]
            try:
                vals = np.array([float(x) for x in row[1:]], dtype=float)
            except ValueError:
                continue
            expr[gene] = vals
    return samples, expr


def rankdata(values):
    values = list(values)
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return np.array(ranks, dtype=float)


def rank_score(samples, expr, genes):
    measured = [g for g in genes if g in expr]
    if not measured:
        return np.full(len(samples), np.nan), []
    all_genes = list(expr)
    matrix = np.vstack([expr[g] for g in all_genes])
    ranks = np.zeros_like(matrix)
    for j in range(matrix.shape[1]):
        ranks[:, j] = rankdata(matrix[:, j])
    idx = [all_genes.index(g) for g in measured]
    raw = ranks[idx, :].mean(axis=0)
    return (raw - 1.0) / max(len(all_genes) - 1.0, 1.0), measured


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def cox_univariate(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    score = np.asarray(score, dtype=float)
    keep = np.isfinite(time) & np.isfinite(score) & (time > 0) & ((event == 0) | (event == 1))
    time, event, score = time[keep], event[keep], score[keep]
    if len(time) < 20 or event.sum() < 5 or score.std(ddof=1) == 0:
        return len(time), int(event.sum()), "", "", "", "", "not_estimable"
    x = (score - score.mean()) / score.std(ddof=1)
    beta = 0.0
    for _ in range(80):
        u = 0.0
        info = 0.0
        for ti, ei, xi in zip(time, event, x):
            if ei != 1:
                continue
            risk = time >= ti
            xr = x[risk]
            w = np.exp(beta * xr)
            sw = w.sum()
            mx = (w * xr).sum() / sw
            mx2 = (w * xr * xr).sum() / sw
            u += xi - mx
            info += mx2 - mx * mx
        if info <= 0:
            break
        step = u / info
        beta += step
        if abs(step) < 1e-8:
            break
    info = 0.0
    for ti, ei in zip(time, event):
        if ei != 1:
            continue
        xr = x[time >= ti]
        w = np.exp(beta * xr)
        sw = w.sum()
        mx = (w * xr).sum() / sw
        mx2 = (w * xr * xr).sum() / sw
        info += mx2 - mx * mx
    se = math.sqrt(1.0 / info) if info > 0 else float("nan")
    p = 2 * norm_sf(abs(beta / se)) if se == se and se > 0 else ""
    return len(time), int(event.sum()), math.exp(beta), math.exp(beta - 1.96 * se), math.exp(beta + 1.96 * se), p, "higher_score_higher_hazard" if beta > 0 else "higher_score_lower_hazard"


def pearson_test(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 5:
        return len(x), "", "", "not_estimable"
    r = np.corrcoef(x, y)[0, 1]
    if abs(r) >= 1:
        p = ""
    else:
        z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(len(x) - 3)
        p = 2 * norm_sf(abs(z))
    return len(x), r, p, "positive" if r > 0 else "negative"


def logistic_numeric(y, score):
    y = np.asarray(y, dtype=float)
    score = np.asarray(score, dtype=float)
    keep = np.isfinite(y) & np.isfinite(score)
    y, score = y[keep], score[keep]
    if len(y) < 20 or y.sum() < 5 or y.sum() > len(y) - 5:
        return len(y), int(y.sum()), "", "", "", "", "not_estimable"
    x = (score - score.mean()) / score.std(ddof=1)
    X = np.column_stack([np.ones(len(x)), x])
    beta = np.zeros(2)
    for _ in range(60):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-8, None)
        info = X.T @ (X * W[:, None])
        grad = X.T @ (y - p)
        step = np.linalg.solve(info, grad)
        beta += step
        if np.max(np.abs(step)) < 1e-8:
            break
    eta = X @ beta
    p = 1 / (1 + np.exp(-eta))
    W = np.clip(p * (1 - p), 1e-8, None)
    info = X.T @ (X * W[:, None])
    cov = np.linalg.inv(info)
    se = math.sqrt(cov[1, 1])
    pval = 2 * norm_sf(abs(beta[1] / se))
    return len(y), int(y.sum()), math.exp(beta[1]), math.exp(beta[1] - 1.96 * se), math.exp(beta[1] + 1.96 * se), pval, "higher_score_higher_odds" if beta[1] > 0 else "higher_score_lower_odds"


def parse_stage_value(x):
    if not x:
        return None
    s = str(x).upper()
    if "STAGE IV" in s or s in {"IV", "IVA", "IVB"}:
        return 4
    if "STAGE III" in s or s in {"III", "IIIA", "IIIB", "IIIC"}:
        return 3
    if "STAGE II" in s or s in {"II", "IIA", "IIB", "IIC"}:
        return 2
    if "STAGE I" in s or s in {"I", "IA", "IB", "IC"}:
        return 1
    return None


def parse_grade_value(x):
    if not x:
        return None
    s = str(x).upper()
    m = re.search(r"G([1-4])", s)
    if m:
        return int(m.group(1))
    return None


def parse_survival_event(x):
    if x in {"1", "0"}:
        return int(x)
    s = str(x).lower()
    if s in {"dead", "deceased", "event"}:
        return 1
    if s in {"alive", "censored"}:
        return 0
    return None


def parse_float(x):
    try:
        if x is None or x == "" or str(x).lower() == "nan":
            return None
        return float(x)
    except ValueError:
        return None


def make_scores(dataset, samples, expr, component_genes, unique_genes):
    score_arrays = {}
    measured = {}
    for comp in COMPONENTS:
        score_arrays[comp], measured[comp] = rank_score(samples, expr, component_genes[comp])
    score_arrays["total_ecotype_unique_gene_score"], measured["total_unique"] = rank_score(samples, expr, unique_genes)
    score_arrays[PRIMARY_SCORE] = np.vstack([score_arrays[c] for c in COMPONENTS]).mean(axis=0)
    return score_arrays, measured


def fmt(x):
    if x == "" or x is None:
        return ""
    if isinstance(x, str):
        return x
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return ""
    return f"{x:.6g}"


def main():
    A05.mkdir(parents=True, exist_ok=True)
    raw_paths = {k: download(v, A05 / (k + (".gz" if v.endswith(".gz") else ".txt"))) for k, v in URLS.items()}

    geneset = read_csv(A02 / "A02_final_ecotype_gene_set.csv")
    component_genes = defaultdict(list)
    gene_components = defaultdict(list)
    for r in geneset:
        component_genes[r["component"]].append(r["gene_symbol"])
        gene_components[r["gene_symbol"]].append(r["component"])
    unique_genes = sorted(gene_components)

    # GSE76427
    gse_meta, gse_header, gse_rows = read_series_matrix(raw_paths["GSE76427_series_matrix"])
    gse_samples = [x.strip('"') for x in gse_header[1:]]
    probe_map = parse_gpl(raw_paths["GPL10558_soft"])
    gse_expr = collapse_probe_rows(gse_rows, gse_samples, probe_map)
    gse_scores, gse_measured = make_scores("GSE76427", gse_samples, gse_expr, component_genes, unique_genes)
    gse_titles = gse_meta.get("!Sample_title", [[]])[0]
    gse_chars = defaultdict(dict)
    for vals in gse_meta.get("!Sample_characteristics_ch1", []):
        for i, v in enumerate(vals):
            if ":" in v and i < len(gse_samples):
                k, val = v.split(":", 1)
                gse_chars[gse_samples[i]][k.strip().lower().replace(" ", "_")] = val.strip()

    # TCGA
    tcga_samples, tcga_expr = load_tcga_hiseq(raw_paths["TCGA_LIHC_HiSeqV2"])
    tcga_scores, tcga_measured = make_scores("TCGA_LIHC", tcga_samples, tcga_expr, component_genes, unique_genes)
    tcga_clin = {}
    with open_text(raw_paths["TCGA_LIHC_clinical"]) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            sid = r.get("sampleID") or r.get("sample") or r.get("Sample")
            if sid:
                tcga_clin[sid] = r
    tcga_surv = {}
    with open_text(raw_paths["TCGA_LIHC_survival"]) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            sid = r.get("sample") or r.get("sampleID") or r.get("Sample")
            if sid:
                tcga_surv[sid] = r

    # Outputs: gene availability
    availability_rows = []
    for dataset, expr, measured in [("GSE76427", gse_expr, gse_measured), ("TCGA_LIHC", tcga_expr, tcga_measured)]:
        for comp in COMPONENTS:
            for g in component_genes[comp]:
                availability_rows.append([dataset, comp, g, "yes" if g in expr else "no"])
        for g in unique_genes:
            availability_rows.append([dataset, "total_ecotype_unique_gene_list", g, "yes" if g in expr else "no"])
    write_csv("A05_gene_availability.csv", ["dataset", "component", "gene_symbol", "available"], availability_rows)

    score_rows = []
    for i, s in enumerate(gse_samples):
        ch = gse_chars.get(s, {})
        score_rows.append([
            "GSE76427", s, gse_titles[i] if i < len(gse_titles) else "", "", "", ch.get("stage", ch.get("tumor_stage", "")),
            ch.get("grade", ch.get("tumor_grade", "")), "", "", "", "",
            gse_scores["malignant_hypoxia_adaptive"][i], gse_scores["endothelial_angiogenesis_interaction"][i],
            gse_scores["spp1_mmp9_tam_matrix_angiogenic_support"][i], gse_scores[PRIMARY_SCORE][i],
            gse_scores["total_ecotype_unique_gene_score"][i],
        ])
    for i, s in enumerate(tcga_samples):
        clin = tcga_clin.get(s, {})
        surv = tcga_surv.get(s, {})
        score_rows.append([
            "TCGA_LIHC", s, "", clin.get("_PATIENT", ""), clin.get("sample_type", ""),
            clin.get("pathologic_stage", clin.get("clinical_stage", "")),
            clin.get("neoplasm_histologic_grade", clin.get("histological_grade", "")),
            clin.get("vascular_tumor_cell_type", clin.get("vascular_invasion", "")),
            surv.get("OS.time", surv.get("OS.time", "")), surv.get("OS", ""),
            surv.get("DSS.time", ""), 
            tcga_scores["malignant_hypoxia_adaptive"][i], tcga_scores["endothelial_angiogenesis_interaction"][i],
            tcga_scores["spp1_mmp9_tam_matrix_angiogenic_support"][i], tcga_scores[PRIMARY_SCORE][i],
            tcga_scores["total_ecotype_unique_gene_score"][i],
        ])
    score_header = [
        "dataset", "sample_id", "sample_title", "patient_id", "sample_type", "stage", "grade", "vascular_annotation",
        "os_time", "os_event", "dss_time", "malignant_hypoxia_adaptive_score",
        "endothelial_angiogenesis_interaction_score", "spp1_mmp9_tam_matrix_angiogenic_support_score",
        PRIMARY_SCORE, "total_ecotype_unique_gene_score",
    ]
    write_csv("A05_ecotype_scores.csv", score_header, score_rows)

    # Disease-context associations: stage/grade/vascular if available.
    disease_rows = []
    for dataset in ["GSE76427", "TCGA_LIHC"]:
        rows = [dict(zip(score_header, r)) for r in score_rows if r[0] == dataset]
        score = [float(r[PRIMARY_SCORE]) for r in rows]
        stage_vals = [parse_stage_value(r["stage"]) for r in rows]
        grade_vals = [parse_grade_value(r["grade"]) for r in rows]
        for field, vals in [("stage_numeric", stage_vals), ("grade_numeric", grade_vals)]:
            n, r, p, direction = pearson_test([v if v is not None else float("nan") for v in vals], score)
            disease_rows.append([dataset, field, PRIMARY_SCORE, n, "", fmt(r), "", "", fmt(p), direction, "secondary_disease_context"])
        high_stage = [1 if v is not None and v >= 3 else (0 if v is not None else float("nan")) for v in stage_vals]
        n, events, orv, lo, hi, p, direction = logistic_numeric(high_stage, score)
        disease_rows.append([dataset, "high_stage_III_IV_vs_I_II", PRIMARY_SCORE, n, events, "", fmt(orv), f"{fmt(lo)} to {fmt(hi)}", fmt(p), direction, "secondary_disease_context"])
    write_csv(
        "A05_disease_context_association_results.csv",
        ["dataset", "disease_context_field", "score_name", "n", "events_or_cases", "correlation", "OR_per_1SD", "OR_95CI", "p_value", "direction", "analysis_scope"],
        disease_rows,
    )

    survival_rows = []
    for dataset in ["GSE76427", "TCGA_LIHC"]:
        rows = [dict(zip(score_header, r)) for r in score_rows if r[0] == dataset]
        if dataset == "TCGA_LIHC":
            time = [parse_float(r["os_time"]) for r in rows]
            event = [parse_survival_event(r["os_event"]) for r in rows]
            score = [float(r[PRIMARY_SCORE]) for r in rows]
            n, ev, hr, lo, hi, p, direction = cox_univariate(
                [t if t is not None else float("nan") for t in time],
                [e if e is not None else -1 for e in event],
                score,
            )
            survival_rows.append([dataset, "OS", PRIMARY_SCORE, n, ev, fmt(hr), fmt(lo), fmt(hi), fmt(p), direction, "secondary_supportive_only"])
        else:
            survival_rows.append([dataset, "OS/RFS", PRIMARY_SCORE, 0, 0, "", "", "", "", "not_available_from_parsed_series_matrix", "field_unavailable"])
    write_csv(
        "A05_survival_or_recurrence_supportive_results.csv",
        ["dataset", "endpoint", "score_name", "n", "events", "HR_per_1SD", "HR_95CI_low", "HR_95CI_high", "p_value", "direction", "analysis_scope"],
        survival_rows,
    )

    robust_rows = []
    for dataset, scores in [("GSE76427", gse_scores), ("TCGA_LIHC", tcga_scores)]:
        for comp in COMPONENTS:
            n, r, p, direction = pearson_test(scores[comp], scores[PRIMARY_SCORE])
            robust_rows.append([dataset, "component_vs_total_correlation", comp, n, fmt(r), fmt(p), direction])
        for leave_out in COMPONENTS:
            kept = [c for c in COMPONENTS if c != leave_out]
            loo = np.vstack([scores[c] for c in kept]).mean(axis=0)
            n, r, p, direction = pearson_test(loo, scores[PRIMARY_SCORE])
            robust_rows.append([dataset, "leave_one_component_out_vs_total", leave_out, n, fmt(r), fmt(p), direction])
        n, r, p, direction = pearson_test(scores["total_ecotype_unique_gene_score"], scores[PRIMARY_SCORE])
        robust_rows.append([dataset, "duplicate_gene_handling_unique_vs_equal_component", "VEGFA_duplicate_once_robustness", n, fmt(r), fmt(p), direction])
    write_csv("A05_component_and_robustness_results.csv", ["dataset", "check_type", "comparison", "n", "estimate", "p_value", "direction"], robust_rows)

    audit = f"""# A05 dataset and endpoint audit

A05 is secondary disease-context externalization only. It does not modify the A02 gene set and does not redefine the project as a pan-HCC prognostic model.

## Datasets downloaded / parsed

- GSE76427 Series Matrix: `{raw_paths['GSE76427_series_matrix'].name}`
- GPL10558 annotation: `{raw_paths['GPL10558_soft'].name}`
- TCGA-LIHC expression: `{raw_paths['TCGA_LIHC_HiSeqV2'].name}`
- TCGA-LIHC clinical matrix: `{raw_paths['TCGA_LIHC_clinical'].name}`
- TCGA-LIHC survival: `{raw_paths['TCGA_LIHC_survival'].name}`

## Sample / gene audit

- GSE76427 samples parsed: {len(gse_samples)}
- GSE76427 gene symbols after probe collapse: {len(gse_expr)}
- TCGA-LIHC samples parsed: {len(tcga_samples)}
- TCGA-LIHC gene symbols parsed: {len(tcga_expr)}

## Endpoint / disease-context fields

- GSE76427: Series Matrix and platform annotation were parsed. Stage/grade/survival fields were not reliably available from the parsed Series Matrix, so GSE76427 is used mainly for score recoverability/gene availability unless fields are manually supplemented later.
- TCGA-LIHC: expression, clinical matrix, and survival table were parsed. Stage/grade fields are used for disease-context association where parseable. OS is used only as secondary supportive evidence.

## Boundary

- No outcome-guided gene pruning.
- No LASSO/RF/SVM/AUC/ROC/cutoff/train-test/nomogram.
- TCGA/GSE76427 do not replace the GSE104580 TACE response endpoint.
"""
    (A05 / "A05_dataset_and_endpoint_audit.md").write_text(audit, encoding="utf-8")

    tcga_os = next(r for r in survival_rows if r[0] == "TCGA_LIHC" and r[1] == "OS")
    verdict = "READY_FOR_A06" if tcga_os[5] != "" or any(r[3] not in {"0", 0} for r in disease_rows) else "BLOCKED_NEED_REPAIR"
    summary = f"""# A05 secondary externalization summary

A05 used the fixed A02 gene set and A03 scoring rule to assess secondary HCC disease-context recoverability. It is not a new prognostic-model phase.

## Chain position

- A02: single-cell/spatial anchor.
- A03: GSE104580 TACE non-response association.
- A04: GSE14520 post-TACE outcome validation with limited specificity.
- A05: secondary disease-context externalization only.

## Interpretation

A05 can support that the ecotype is recoverable in broader HCC bulk data and may reflect adverse tumor biology. It must not become the main OS/RFS story.

Final verdict: `{verdict}`.
"""
    (A05 / "A05_secondary_externalization_summary.md").write_text(summary, encoding="utf-8")
    final = f"""# final_A05_verdict

`{verdict}`

## Boundary

A05 is secondary disease-context externalization only. It does not alter A02/A03/A04 conclusions and must not reframe the project as a generic HCC prognostic model.

## A06-ready files

- Scores: `results/a05_secondary_externalization/A05_ecotype_scores.csv`
- Disease-context associations: `results/a05_secondary_externalization/A05_disease_context_association_results.csv`
- Supportive survival/recurrence: `results/a05_secondary_externalization/A05_survival_or_recurrence_supportive_results.csv`

## Still forbidden

- No edits to A02 gene set.
- No outcome-guided gene pruning.
- No LASSO/RF/SVM/AUC/ROC/cutoff/train-test/nomogram.
- No OS/RFS-led reframing.
- No TCGA/GSE76427 dominance over the endpoint-first TACE response chain.
- No ICI/TLS/MVI/perivascular/immune-exclusion/early-recurrence pivot.
"""
    (A05 / "final_A05_verdict.md").write_text(final, encoding="utf-8")


if __name__ == "__main__":
    main()
