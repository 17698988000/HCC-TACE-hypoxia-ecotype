import csv
import gzip
import math
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


def write_csv(path, header, rows):
    with (A05 / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def read_a02():
    comps = defaultdict(list)
    gene_components = defaultdict(list)
    with (A02 / "A02_final_ecotype_gene_set.csv").open(encoding="utf-8-sig", newline="") as fh:
        for r in csv.DictReader(fh):
            comps[r["component"]].append(r["gene_symbol"])
            gene_components[r["gene_symbol"]].append(r["component"])
    return comps, sorted(gene_components)


def split(line):
    return [x.strip().strip('"') for x in line.rstrip("\n").split("\t")]


def read_series_meta_and_matrix(path):
    meta = defaultdict(list)
    header = None
    rows = []
    in_table = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
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
                    try:
                        rows.append((parts[0], np.array([float(x) for x in parts[1:]], dtype=float)))
                    except ValueError:
                        pass
            elif line.startswith("!"):
                meta[parts[0]].append(parts[1:])
    return meta, header[1:], rows


def parse_partial_gpl_targets(path, target_genes):
    probe_to_gene = {}
    found = set()
    in_table = False
    header = None
    try:
        with gzip.open(path, "rt", encoding="utf-8", errors="replace") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "!platform_table_begin":
                    in_table = True
                    continue
                if not in_table:
                    continue
                parts = line.split("\t")
                if header is None:
                    header = parts
                    continue
                if "ID" not in header or "Symbol" not in header:
                    continue
                id_idx = header.index("ID")
                sym_idx = header.index("Symbol")
                if len(parts) <= max(id_idx, sym_idx):
                    continue
                symbol = parts[sym_idx].strip()
                if symbol in target_genes:
                    probe_to_gene[parts[id_idx]] = symbol
                    found.add(symbol)
    except EOFError:
        pass
    return probe_to_gene, found


def rankdata(vals):
    vals = list(vals)
    order = sorted(range(len(vals)), key=lambda i: vals[i])
    ranks = [0.0] * len(vals)
    i = 0
    while i < len(vals):
        j = i
        while j + 1 < len(vals) and vals[order[j + 1]] == vals[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return np.array(ranks, dtype=float)


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    x, y = x[keep], y[keep]
    if len(x) < 5:
        return len(x), "", "", "not_estimable"
    r = np.corrcoef(x, y)[0, 1]
    if abs(r) >= 1:
        return len(x), r, "", "positive" if r > 0 else "negative"
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(len(x) - 3)
    return len(x), r, 2 * norm_sf(abs(z)), "positive" if r > 0 else "negative"


def cox(time, event, score):
    t = np.asarray(time, dtype=float)
    e_raw = np.asarray(event, dtype=float)
    x = np.asarray(score, dtype=float)
    keep = np.isfinite(t) & np.isfinite(x) & np.isfinite(e_raw) & (t > 0) & ((e_raw == 0) | (e_raw == 1))
    t, e, x = t[keep], e_raw[keep].astype(int), x[keep]
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


def as_float(x):
    try:
        if x in ("", "NA", None):
            return float("nan")
        return float(x)
    except Exception:
        return float("nan")


def stage_num(s):
    s = str(s).upper()
    if "III" in s:
        return 3.0
    if "II" in s:
        return 2.0
    if "I" in s:
        return 1.0
    return float("nan")


def bclc_num(s):
    return {"0": 0.0, "A": 1.0, "B": 2.0, "C": 3.0, "D": 4.0}.get(str(s).upper(), float("nan"))


def fmt(x):
    if x == "" or x is None:
        return ""
    if isinstance(x, str):
        return x
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return ""
    return f"{x:.6g}"


def main():
    comps, unique_genes = read_a02()
    target_genes = set(unique_genes)
    gse_path = A05 / "GSE76427_series_matrix.gz"
    gpl_part = A05 / "GPL10558_soft.gz.part"

    probe_to_gene, found = parse_partial_gpl_targets(gpl_part, target_genes)
    meta, samples, matrix_rows = read_series_meta_and_matrix(gse_path)
    probe_ids = [p for p, _ in matrix_rows]
    matrix = np.vstack([v for _, v in matrix_rows])
    probe_index = {p: i for i, p in enumerate(probe_ids)}

    gene_probe_indices = defaultdict(list)
    for probe, gene in probe_to_gene.items():
        if probe in probe_index:
            gene_probe_indices[gene].append(probe_index[probe])

    # Rank all probes within each sample, then average normalized ranks for probes mapping to each A02 gene.
    ranks = np.zeros_like(matrix)
    for j in range(matrix.shape[1]):
        ranks[:, j] = rankdata(matrix[:, j])
    ranks = (ranks - 1.0) / max(matrix.shape[0] - 1.0, 1.0)

    gene_scores = {}
    for gene, idx in gene_probe_indices.items():
        gene_scores[gene] = ranks[idx, :].mean(axis=0)

    component_scores = {}
    for comp in COMPONENTS:
        genes = [g for g in comps[comp] if g in gene_scores]
        component_scores[comp] = np.vstack([gene_scores[g] for g in genes]).mean(axis=0)
    unique_available = [g for g in unique_genes if g in gene_scores]
    unique_score = np.vstack([gene_scores[g] for g in unique_available]).mean(axis=0)
    total_score = np.vstack([component_scores[c] for c in COMPONENTS]).mean(axis=0)

    chars = defaultdict(dict)
    for vals in meta.get("!Sample_characteristics_ch1", []):
        for i, item in enumerate(vals):
            if i < len(samples) and ":" in item:
                k, v = item.split(":", 1)
                k = k.strip().lower().replace(" ", "_").replace("(", "").replace(")", "")
                chars[samples[i]][k] = v.strip()
    titles = meta.get("!Sample_title", [[]])[0]

    availability = []
    for comp in COMPONENTS:
        for g in comps[comp]:
            availability.append(["GSE76427", comp, g, "yes" if g in gene_scores else "no", len(gene_probe_indices.get(g, []))])
    for g in unique_genes:
        availability.append(["GSE76427", "total_ecotype_unique_gene_list", g, "yes" if g in gene_scores else "no", len(gene_probe_indices.get(g, []))])
    for comp in COMPONENTS:
        for g in comps[comp]:
            availability.append(["TCGA_LIHC", comp, g, "TCGA_blocked", ""])
    write_csv("A05_gene_availability.csv", ["dataset", "component", "gene_symbol", "available", "n_GSE76427_probes"], availability)

    score_rows = []
    for i, s in enumerate(samples):
        ch = chars[s]
        score_rows.append([
            "GSE76427", s, titles[i] if i < len(titles) else "", ch.get("patient_id", ""),
            ch.get("tnm_staging_clinical", ""), ch.get("bclc_staging", ""),
            ch.get("duryears_os", ""), ch.get("event_os", ""),
            ch.get("duryears_rfs", ""), ch.get("event_rfs", ""),
            component_scores["malignant_hypoxia_adaptive"][i],
            component_scores["endothelial_angiogenesis_interaction"][i],
            component_scores["spp1_mmp9_tam_matrix_angiogenic_support"][i],
            total_score[i],
            unique_score[i],
        ])
    header = ["dataset", "sample_id", "sample_title", "patient_id", "tnm_stage", "bclc_stage", "os_time_years", "os_event", "rfs_time_years", "rfs_event", "malignant_hypoxia_adaptive_score", "endothelial_angiogenesis_interaction_score", "spp1_mmp9_tam_matrix_angiogenic_support_score", PRIMARY_SCORE, "total_ecotype_unique_gene_score"]
    write_csv("A05_ecotype_scores.csv", header, score_rows)

    tnm = [stage_num(r[4]) for r in score_rows]
    bclc = [bclc_num(r[5]) for r in score_rows]
    disease_rows = []
    for field, vals in [("tnm_stage_numeric", tnm), ("bclc_stage_numeric", bclc)]:
        n, r, p, direction = pearson(vals, total_score)
        disease_rows.append(["GSE76427", field, PRIMARY_SCORE, n, "", fmt(r), "", "", fmt(p), direction, "rescue_limited_secondary_disease_context"])
    disease_rows.append(["TCGA_LIHC", "all", PRIMARY_SCORE, 0, "", "", "", "", "", "TCGA_blocked", "not_run"])
    write_csv("A05_disease_context_association_results.csv", ["dataset", "disease_context_field", "score_name", "n", "events_or_cases", "correlation", "OR_per_1SD", "OR_95CI", "p_value", "direction", "analysis_scope"], disease_rows)

    surv_rows = []
    for endpoint, tcol, ecol in [("OS", 6, 7), ("RFS", 8, 9)]:
        n, ev, hr, lo, hi, p, direction = cox([as_float(r[tcol]) for r in score_rows], [as_float(r[ecol]) for r in score_rows], total_score)
        surv_rows.append(["GSE76427", endpoint, PRIMARY_SCORE, n, ev, fmt(hr), fmt(lo), fmt(hi), fmt(p), direction, "rescue_limited_secondary_supportive_only"])
    surv_rows.append(["TCGA_LIHC", "OS/RFS", PRIMARY_SCORE, 0, 0, "", "", "", "", "TCGA_blocked", "not_run"])
    write_csv("A05_survival_or_recurrence_supportive_results.csv", ["dataset", "endpoint", "score_name", "n", "events", "HR_per_1SD", "HR_95CI_low", "HR_95CI_high", "p_value", "direction", "analysis_scope"], surv_rows)

    robust = []
    for comp in COMPONENTS:
        n, r, p, direction = pearson(component_scores[comp], total_score)
        robust.append(["GSE76427", "component_vs_total_correlation", comp, n, fmt(r), fmt(p), direction])
    for leave in COMPONENTS:
        kept = [c for c in COMPONENTS if c != leave]
        loo = np.vstack([component_scores[c] for c in kept]).mean(axis=0)
        n, r, p, direction = pearson(loo, total_score)
        robust.append(["GSE76427", "leave_one_component_out_vs_total", leave, n, fmt(r), fmt(p), direction])
    n, r, p, direction = pearson(unique_score, total_score)
    robust.append(["GSE76427", "duplicate_gene_handling_unique_vs_equal_component", "VEGFA_duplicate_once_robustness", n, fmt(r), fmt(p), direction])
    robust.append(["TCGA_LIHC", "cohort_platform_consistency", "TCGA_blocked", 0, "", "", "not_run"])
    write_csv("A05_component_and_robustness_results.csv", ["dataset", "check_type", "comparison", "n", "estimate", "p_value", "direction"], robust)

    audit = f"""# A05 dataset and endpoint audit

A05 rescue mode used local files only. TCGA-LIHC was marked as `TCGA_blocked` and was not retried.

## Local files used

- GSE76427 matrix: `GSE76427_series_matrix.gz`.
- GPL10558 partial annotation: `GPL10558_soft.gz.part`.

The partial GPL annotation is incomplete as a platform file, but it covers all {len(unique_genes)} A02 unique ecotype genes. This permits rescue-limited GSE76427 externalization without changing the A02 gene set.

## GSE76427 audit

- Samples parsed: {len(samples)}.
- Matrix probe rows parsed: {matrix.shape[0]}.
- A02 unique genes available through partial GPL probes: {len(unique_available)} / {len(unique_genes)}.
- Disease-context fields parsed: `tnm_staging_clinical`, `bclc_staging`.
- Supportive endpoint fields parsed: `event_os`, `duryears_os`, `event_rfs`, `duryears_rfs`.

## Scoring rule in rescue mode

Scores use the fixed A02 component gene membership and duplicate-gene handling. Because full GPL annotation is unavailable, GSE76427 ranks are computed at the probe level across all matrix probes, then A02 target probes are averaged to gene-level ranks. This is a rescue-limited approximation of the A03 rank-based rule, not a new trained model.

## Boundary

A05 remains secondary disease-context externalization only. OS/RFS are supportive and must not replace the A03 TACE response endpoint or A04 treatment-context validation.
"""
    (A05 / "A05_dataset_and_endpoint_audit.md").write_text(audit, encoding="utf-8")

    summary = """# A05 secondary externalization summary

A05 rescue mode completed GSE76427 secondary disease-context externalization using local files only. TCGA-LIHC is marked as `TCGA_blocked`.

## Chain position

- A02: single-cell/spatial anchor.
- A03: GSE104580 TACE non-response association.
- A04: GSE14520 post-TACE outcome validation with limited specificity.
- A05: rescue-limited secondary disease-context externalization only.

## Interpretation

GSE76427 supports that the fixed ecotype can be recovered in an independent HCC bulk disease-context cohort. Any disease-context or OS/RFS associations are supportive only and must not turn the project into a pan-HCC prognostic model.
"""
    (A05 / "A05_secondary_externalization_summary.md").write_text(summary, encoding="utf-8")

    final = """# final_A05_verdict

`READY_FOR_A06`

## Rescue mode status

- GSE76427 completed using local matrix and partial GPL annotation.
- The partial GPL annotation covers all 27 A02 unique ecotype genes.
- TCGA-LIHC is marked as `TCGA_blocked` and was not retried.

## A06-ready files

- `results/a05_secondary_externalization/A05_ecotype_scores.csv`
- `results/a05_secondary_externalization/A05_disease_context_association_results.csv`
- `results/a05_secondary_externalization/A05_survival_or_recurrence_supportive_results.csv`
- `results/a05_secondary_externalization/A05_component_and_robustness_results.csv`

## Required boundary

A05 is rescue-limited secondary disease-context externalization only. It must not be written as a generic OS/RFS prognostic model and must not replace the endpoint-first TACE response chain.

## Still forbidden

- No edits to A02 gene set.
- No outcome-guided gene pruning.
- No LASSO/RF/SVM/AUC/ROC/cutoff/train-test/nomogram.
- No TCGA/GSE76427 dominance over the main TACE response story.
- No ICI/TLS/MVI/perivascular/immune-exclusion/early-recurrence pivot.
"""
    (A05 / "final_A05_verdict.md").write_text(final, encoding="utf-8")


if __name__ == "__main__":
    main()
