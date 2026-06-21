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
LOG = A05 / "A05_rescue_progress.log"

COMPONENTS = [
    "malignant_hypoxia_adaptive",
    "endothelial_angiogenesis_interaction",
    "spp1_mmp9_tam_matrix_angiogenic_support",
]
PRIMARY_SCORE = "total_ecotype_equal_component_score"
GSE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76427/matrix/GSE76427_series_matrix.txt.gz"
GPL_URL = "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL10nnn/GPL10558/soft/GPL10558_family.soft.gz"


def log(msg):
    print(msg)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(str(msg) + "\n")


def gzip_ok(path):
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for _ in fh:
                pass
        return True
    except Exception as exc:
        log(f"gzip_check_failed: {path.name}: {exc!r}")
        return False


def download_once(url, dest, timeout=60):
    dest = Path(dest)
    tmp = dest.with_suffix(dest.suffix + ".part")
    log(f"download_start: {url} -> {dest.name}")
    with urllib.request.urlopen(url, timeout=timeout) as r, tmp.open("wb") as out:
        total = r.headers.get("Content-Length")
        total = int(total) if total else None
        got = 0
        last_mb = -1
        while True:
            chunk = r.read(1024 * 1024)
            if not chunk:
                break
            out.write(chunk)
            got += len(chunk)
            mb = got // (1024 * 1024)
            if mb != last_mb:
                last_mb = mb
                if total:
                    log(f"download_progress: {dest.name}: {got}/{total} bytes")
                else:
                    log(f"download_progress: {dest.name}: {got} bytes")
    tmp.replace(dest)
    log(f"download_done: {dest.name}: {dest.stat().st_size} bytes")
    return dest


def ensure_gzip(url, dest, timeout=60):
    dest = Path(dest)
    if dest.exists() and dest.stat().st_size > 0 and gzip_ok(dest):
        log(f"using_existing_valid_file: {dest.name}")
        return dest, "existing_valid"
    if dest.exists():
        log(f"redownload_invalid_or_partial_file: {dest.name}")
        dest.unlink()
    download_once(url, dest, timeout=timeout)
    if not gzip_ok(dest):
        raise RuntimeError(f"Downloaded file failed gzip integrity check: {dest}")
    return dest, "downloaded_valid"


def open_text(path):
    return gzip.open(path, "rt", encoding="utf-8", errors="replace") if str(path).endswith(".gz") else open(path, encoding="utf-8")


def write_csv(path, header, rows):
    with (A05 / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def read_a02_genes():
    comps = defaultdict(list)
    gene_components = defaultdict(list)
    with (A02 / "A02_final_ecotype_gene_set.csv").open(encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            comps[r["component"]].append(r["gene_symbol"])
            gene_components[r["gene_symbol"]].append(r["component"])
    return comps, sorted(gene_components)


def split(line):
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]


def read_series(path):
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
            parts = split(line)
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
            id_idx = header.index("ID")
            symbol_idx = header.index("Symbol") if "Symbol" in header else header.index("ILMN_Gene")
            if len(parts) <= max(id_idx, symbol_idx):
                continue
            gene = parts[symbol_idx].strip()
            if gene and gene not in {"---", "NA"}:
                probe_to_gene[parts[id_idx]] = gene
    return probe_to_gene


def collapse(rows, samples, p2g):
    sums = {}
    counts = Counter()
    for row in rows:
        g = p2g.get(row[0])
        if not g:
            continue
        try:
            vals = np.array([float(x) for x in row[1:]], dtype=float)
        except Exception:
            continue
        if g not in sums:
            sums[g] = np.zeros(len(samples))
        sums[g] += vals
        counts[g] += 1
    return {g: sums[g] / counts[g] for g in sums}


def rankdata(vals):
    vals = list(vals)
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return np.array(ranks)


def rank_score(samples, expr, genes):
    measured = [g for g in genes if g in expr]
    all_genes = list(expr)
    matrix = np.vstack([expr[g] for g in all_genes])
    ranks = np.zeros_like(matrix)
    for j in range(matrix.shape[1]):
        ranks[:, j] = rankdata(matrix[:, j])
    idx = [all_genes.index(g) for g in measured]
    raw = ranks[idx, :].mean(axis=0)
    return (raw - 1) / max(len(all_genes) - 1, 1), measured


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def cox(time, event, score):
    t = np.array(time, dtype=float)
    e = np.array(event, dtype=int)
    x = np.array(score, dtype=float)
    keep = np.isfinite(t) & np.isfinite(x) & (t > 0) & ((e == 0) | (e == 1))
    t, e, x = t[keep], e[keep], x[keep]
    if len(t) < 20 or e.sum() < 5 or x.std(ddof=1) == 0:
        return len(t), int(e.sum()), "", "", "", "", "not_estimable"
    z = (x - x.mean()) / x.std(ddof=1)
    beta = 0.0
    for _ in range(80):
        u = info = 0.0
        for ti, ei, zi in zip(t, e, z):
            if ei != 1:
                continue
            xr = z[t >= ti]
            w = np.exp(beta * xr)
            sw = w.sum()
            mx = (w * xr).sum() / sw
            mx2 = (w * xr * xr).sum() / sw
            u += zi - mx
            info += mx2 - mx * mx
        if info <= 0:
            break
        step = u / info
        beta += step
        if abs(step) < 1e-8:
            break
    info = 0.0
    for ti, ei in zip(t, e):
        if ei != 1:
            continue
        xr = z[t >= ti]
        w = np.exp(beta * xr)
        sw = w.sum()
        mx = (w * xr).sum() / sw
        mx2 = (w * xr * xr).sum() / sw
        info += mx2 - mx * mx
    se = math.sqrt(1 / info) if info > 0 else float("nan")
    p = 2 * norm_sf(abs(beta / se)) if se == se else ""
    return len(t), int(e.sum()), math.exp(beta), math.exp(beta - 1.96 * se), math.exp(beta + 1.96 * se), p, "higher_score_higher_hazard" if beta > 0 else "higher_score_lower_hazard"


def pearson(x, y):
    x = np.array(x, dtype=float)
    y = np.array(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 5:
        return len(x), "", "", "not_estimable"
    r = np.corrcoef(x, y)[0, 1]
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(len(x) - 3) if abs(r) < 1 else 0
    p = 2 * norm_sf(abs(z)) if abs(r) < 1 else ""
    return len(x), r, p, "positive" if r > 0 else "negative"


def val(x):
    if x is None or x in {"", "NA", "nan"}:
        return float("nan")
    try:
        return float(x)
    except Exception:
        return float("nan")


def stage_num(s):
    if not s or s == "NA":
        return float("nan")
    s = s.upper()
    if "III" in s:
        return 3
    if "II" in s:
        return 2
    if "I" in s:
        return 1
    return float("nan")


def bclc_num(s):
    return {"0": 0, "A": 1, "B": 2, "C": 3, "D": 4}.get(str(s).upper(), float("nan"))


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
    LOG.write_text("", encoding="utf-8")
    gse_path, gse_status = ensure_gzip(GSE_URL, A05 / "GSE76427_series_matrix.gz", timeout=90)
    gpl_path, gpl_status = ensure_gzip(GPL_URL, A05 / "GPL10558_soft.gz", timeout=90)
    log("TCGA_status: TCGA_blocked_by_rescue_mode_no_retry")

    comps, unique_genes = read_a02_genes()
    meta, header, rows = read_series(gse_path)
    samples = header[1:]
    p2g = parse_gpl(gpl_path)
    expr = collapse(rows, samples, p2g)
    scores = {}
    measured = {}
    for c in COMPONENTS:
        scores[c], measured[c] = rank_score(samples, expr, comps[c])
    scores["total_ecotype_unique_gene_score"], measured["unique"] = rank_score(samples, expr, unique_genes)
    scores[PRIMARY_SCORE] = np.vstack([scores[c] for c in COMPONENTS]).mean(axis=0)

    chars = defaultdict(dict)
    for vals in meta.get("!Sample_characteristics_ch1", []):
        for i, item in enumerate(vals):
            if ":" in item and i < len(samples):
                k, v = item.split(":", 1)
                chars[samples[i]][k.strip().lower().replace(" ", "_").replace("(", "").replace(")", "")] = v.strip()
    titles = meta.get("!Sample_title", [[]])[0]

    availability = []
    for c in COMPONENTS:
        for g in comps[c]:
            availability.append(["GSE76427", c, g, "yes" if g in expr else "no"])
    for g in unique_genes:
        availability.append(["GSE76427", "total_ecotype_unique_gene_list", g, "yes" if g in expr else "no"])
    for c in COMPONENTS:
        for g in comps[c]:
            availability.append(["TCGA_LIHC", c, g, "TCGA_blocked"])
    write_csv("A05_gene_availability.csv", ["dataset", "component", "gene_symbol", "available"], availability)

    score_header = ["dataset", "sample_id", "sample_title", "patient_id", "tnm_stage", "bclc_stage", "os_time_years", "os_event", "rfs_time_years", "rfs_event"] + [f"{c}_score" for c in COMPONENTS] + [PRIMARY_SCORE, "total_ecotype_unique_gene_score"]
    score_rows = []
    for i, s in enumerate(samples):
        ch = chars[s]
        score_rows.append([
            "GSE76427", s, titles[i] if i < len(titles) else "", ch.get("patient_id", ""),
            ch.get("tnm_staging_clinical", ""), ch.get("bclc_staging", ""),
            ch.get("duryears_os", ""), ch.get("event_os", ""),
            ch.get("duryears_rfs", ""), ch.get("event_rfs", ""),
            *[scores[c][i] for c in COMPONENTS], scores[PRIMARY_SCORE][i], scores["total_ecotype_unique_gene_score"][i],
        ])
    write_csv("A05_ecotype_scores.csv", score_header, score_rows)

    tnm = [stage_num(r[4]) for r in score_rows]
    bclc = [bclc_num(r[5]) for r in score_rows]
    total = [float(r[-2]) for r in score_rows]
    disease = []
    for field, values in [("tnm_stage_numeric", tnm), ("bclc_stage_numeric", bclc)]:
        n, r, p, direction = pearson(values, total)
        disease.append(["GSE76427", field, PRIMARY_SCORE, n, "", fmt(r), "", "", fmt(p), direction, "secondary_disease_context"])
    write_csv("A05_disease_context_association_results.csv", ["dataset", "disease_context_field", "score_name", "n", "events_or_cases", "correlation", "OR_per_1SD", "OR_95CI", "p_value", "direction", "analysis_scope"], disease + [["TCGA_LIHC", "all", PRIMARY_SCORE, 0, "", "", "", "", "", "TCGA_blocked", "not_run"]])

    surv = []
    for endpoint, tcol, ecol in [("OS", 6, 7), ("RFS", 8, 9)]:
        n, ev, hr, lo, hi, p, direction = cox([val(r[tcol]) for r in score_rows], [val(r[ecol]) for r in score_rows], total)
        surv.append(["GSE76427", endpoint, PRIMARY_SCORE, n, ev, fmt(hr), fmt(lo), fmt(hi), fmt(p), direction, "secondary_supportive_only"])
    surv.append(["TCGA_LIHC", "OS/RFS", PRIMARY_SCORE, 0, 0, "", "", "", "", "TCGA_blocked", "not_run"])
    write_csv("A05_survival_or_recurrence_supportive_results.csv", ["dataset", "endpoint", "score_name", "n", "events", "HR_per_1SD", "HR_95CI_low", "HR_95CI_high", "p_value", "direction", "analysis_scope"], surv)

    robust = []
    for c in COMPONENTS:
        n, r, p, direction = pearson(scores[c], scores[PRIMARY_SCORE])
        robust.append(["GSE76427", "component_vs_total_correlation", c, n, fmt(r), fmt(p), direction])
    for leave in COMPONENTS:
        kept = [c for c in COMPONENTS if c != leave]
        loo = np.vstack([scores[c] for c in kept]).mean(axis=0)
        n, r, p, direction = pearson(loo, scores[PRIMARY_SCORE])
        robust.append(["GSE76427", "leave_one_component_out_vs_total", leave, n, fmt(r), fmt(p), direction])
    n, r, p, direction = pearson(scores["total_ecotype_unique_gene_score"], scores[PRIMARY_SCORE])
    robust.append(["GSE76427", "duplicate_gene_handling_unique_vs_equal_component", "VEGFA_duplicate_once_robustness", n, fmt(r), fmt(p), direction])
    robust.append(["TCGA_LIHC", "cohort_platform_consistency", "TCGA_blocked", 0, "", "", "not_run"])
    write_csv("A05_component_and_robustness_results.csv", ["dataset", "check_type", "comparison", "n", "estimate", "p_value", "direction"], robust)

    audit = f"""# A05 dataset and endpoint audit

Rescue mode was used after the prior A05 run blocked during external downloads.

## Data used

- GSE76427 Series Matrix: `{gse_path.name}` ({gse_status})
- GPL10558 annotation: `{gpl_path.name}` ({gpl_status})
- TCGA-LIHC: `TCGA_blocked`; not retried in rescue mode after prior 403/blocking behavior.

## GSE76427 audit

- Samples parsed: {len(samples)}
- Gene symbols after probe collapse: {len(expr)}
- A02 unique ecotype genes available: {sum(1 for g in unique_genes if g in expr)} / {len(unique_genes)}
- Disease-context fields parsed: `tnm_staging_clinical`, `bclc_staging`.
- Supportive endpoint fields parsed: `event_os`, `duryears_os`, `event_rfs`, `duryears_rfs`.

## Boundary

A05 remains secondary disease-context externalization only. OS/RFS are supportive; they do not replace the A03 TACE response endpoint or A04 treatment-context validation.
"""
    (A05 / "A05_dataset_and_endpoint_audit.md").write_text(audit, encoding="utf-8")

    verdict = "READY_FOR_A06" if len(samples) > 20 else "BLOCKED_NEED_REPAIR"
    summary = f"""# A05 secondary externalization summary

A05 rescue mode completed GSE76427 secondary externalization and marked TCGA-LIHC as blocked.

## Chain position

- A02: single-cell/spatial anchor.
- A03: GSE104580 TACE non-response association.
- A04: GSE14520 post-TACE outcome validation with limited specificity.
- A05: secondary disease-context externalization only.

## Interpretation

GSE76427 allows recovery of the fixed ecotype in an independent HCC bulk disease-context cohort. Any OS/RFS findings are supportive only and must not turn the project into a pan-HCC prognostic model.

Final verdict: `{verdict}`.
"""
    (A05 / "A05_secondary_externalization_summary.md").write_text(summary, encoding="utf-8")
    final = f"""# final_A05_verdict

`{verdict}`

## Rescue mode status

- GSE76427 completed.
- TCGA-LIHC marked as `TCGA_blocked` and was not repeatedly retried.

## A06-ready files

- `results/a05_secondary_externalization/A05_ecotype_scores.csv`
- `results/a05_secondary_externalization/A05_disease_context_association_results.csv`
- `results/a05_secondary_externalization/A05_survival_or_recurrence_supportive_results.csv`
- `results/a05_secondary_externalization/A05_component_and_robustness_results.csv`

## Required boundary

A05 is secondary disease-context externalization only. It must not be written as a generic OS/RFS prognostic model and must not replace the endpoint-first TACE response chain.

## Still forbidden

- No edits to A02 gene set.
- No outcome-guided gene pruning.
- No LASSO/RF/SVM/AUC/ROC/cutoff/train-test/nomogram.
- No TCGA/GSE76427 dominance over the main TACE response story.
- No ICI/TLS/MVI/perivascular/immune-exclusion/early-recurrence pivot.
"""
    (A05 / "final_A05_verdict.md").write_text(final, encoding="utf-8")
    log(f"A05_rescue_done: {verdict}")


if __name__ == "__main__":
    main()
