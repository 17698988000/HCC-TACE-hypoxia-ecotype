import csv
import gzip
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
A01 = ROOT / "results" / "a01_data_preparation"
A02 = ROOT / "results" / "a02_ecotype_construction"


PREDEFINED_PROGRAMS = {
    "malignant_hypoxia_adaptive": [
        "CA9", "VEGFA", "SLC2A1", "HK2", "LDHA", "PDK1", "BNIP3", "EGLN3",
        "ADM", "ENO1", "NDRG1", "HILPDA", "DDIT4", "ALDOA",
    ],
    "endothelial_angiogenesis_interaction": [
        "PECAM1", "VWF", "KDR", "FLT1", "ENG", "CD34", "EMCN", "RAMP2",
        "ACKR1", "ESAM", "ANGPT2", "TEK", "PLVAP", "VEGFB", "PGF",
        "LGALS9", "CXCL12", "CXCR4", "DLL4", "NOTCH1", "JAG1",
    ],
    "spp1_mmp9_tam_matrix_angiogenic_support": [
        "SPP1", "MMP9", "CD68", "LST1", "C1QA", "C1QB", "C1QC", "APOE",
        "MARCO", "TREM2", "LGALS3", "MSR1", "VEGFA",
    ],
}


COMPONENT_TARGETS = {
    "GSE149614": {
        "malignant_hypoxia_adaptive": {"celltypes": {"Hepatocyte"}, "sites": {"Tumor", "PVTT"}},
        "endothelial_angiogenesis_interaction": {"celltypes": {"Endothelial"}, "sites": None},
        "spp1_mmp9_tam_matrix_angiogenic_support": {"celltypes": {"Myeloid"}, "sites": None},
    },
    "GSE151530": {
        "malignant_hypoxia_adaptive": {"celltypes": {"Malignant cells"}, "sites": None},
        "endothelial_angiogenesis_interaction": {"celltypes": {"TECs"}, "sites": None},
        "spp1_mmp9_tam_matrix_angiogenic_support": {"celltypes": {"TAMs"}, "sites": None},
    },
}


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_csv(path, header, rows):
    with (A02 / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def load_matrix_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        cols = header[1:]
        data = {}
        for row in reader:
            if not row:
                continue
            data[row[0]] = [float(x) for x in row[1:]]
    return cols, data


def load_matrix_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        cols = header[1:]
        data = {}
        for row in reader:
            if not row:
                continue
            data[row[0]] = [float(x) for x in row[1:]]
    return cols, data


def load_gse277104_qnorm(path):
    with Path(path).open(encoding="utf-8-sig") as fh:
        reader = csv.reader(fh, delimiter="\t")
        cols = next(reader)
        data = {}
        for row in reader:
            if not row:
                continue
            data[row[0]] = [float(x) for x in row[1:]]
    return cols, data


def read_series_matrix_meta(path):
    meta = defaultdict(list)
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line == "!series_matrix_table_begin":
                break
            if line.startswith("!"):
                parts = [x.strip('"') for x in line.split("\t")]
                meta[parts[0]].append(parts[1:])
    return meta


def dcc_name_from_supplementary(url):
    name = url.rsplit("/", 1)[-1].replace(".gz", "")
    if "_" in name and name.startswith("GSM"):
        return name.split("_", 1)[1]
    return name


def mean(vals):
    vals = list(vals)
    return sum(vals) / len(vals) if vals else float("nan")


def median(vals):
    vals = sorted(vals)
    if not vals:
        return float("nan")
    mid = len(vals) // 2
    if len(vals) % 2:
        return vals[mid]
    return (vals[mid - 1] + vals[mid]) / 2


def safe_fmt(x):
    if isinstance(x, str):
        return x
    if x != x:
        return ""
    return f"{x:.4f}"


def parse_group_metadata(dataset):
    if dataset == "GSE149614":
        rows = read_tsv(A01 / "metadata_GSE149614_pseudobulk_group_metadata.tsv")
        out = {}
        for r in rows:
            out[r["group_id"]] = {
                "sample": r["sample"],
                "site": r["site"],
                "celltype": r["celltype"],
                "n_cells": int(r["n_cells"]),
            }
        return out
    rows = read_tsv(A01 / "metadata_GSE151530_HCC_pseudobulk_group_metadata.tsv")
    out = {}
    for r in rows:
        out[r["group_id"]] = {
            "sample": r["S_ID"],
            "site": "HCC",
            "celltype": r["celltype"],
            "n_cells": int(r["n_cells"]),
        }
    return out


def single_cell_anchor(dataset, matrix_file):
    cols, data = load_matrix_csv(A01 / matrix_file)
    meta = parse_group_metadata(dataset)
    per_cell_log = {}
    for gene, vals in data.items():
        converted = []
        for col, val in zip(cols, vals):
            n = max(meta[col]["n_cells"], 1)
            converted.append(math.log1p(val / n))
        per_cell_log[gene] = converted

    rows = []
    component_gene_stats = {}
    for component, genes in PREDEFINED_PROGRAMS.items():
        rules = COMPONENT_TARGETS[dataset][component]
        target_idx = []
        other_idx = []
        for i, col in enumerate(cols):
            m = meta[col]
            is_target = m["celltype"] in rules["celltypes"]
            if rules["sites"] is not None:
                is_target = is_target and m["site"] in rules["sites"]
            if is_target:
                target_idx.append(i)
            else:
                other_idx.append(i)
        for gene in genes:
            if gene not in per_cell_log:
                rows.append([
                    dataset, component, gene, "no", len(target_idx), "", "", "", "not_detected_in_candidate_pseudobulk",
                ])
                continue
            vals = per_cell_log[gene]
            target_vals = [vals[i] for i in target_idx]
            other_vals = [vals[i] for i in other_idx]
            target_mean = mean(target_vals)
            other_mean = mean(other_vals)
            delta = target_mean - other_mean
            target_detected_fraction = sum(1 for i in target_idx if vals[i] > 0) / max(len(target_idx), 1)
            support = "target_enriched" if delta > 0 and target_detected_fraction >= 0.2 else "detected_not_target_enriched"
            component_gene_stats[(dataset, component, gene)] = {
                "target_mean": target_mean,
                "other_mean": other_mean,
                "delta": delta,
                "detected_fraction": target_detected_fraction,
                "support": support,
            }
            rows.append([
                dataset,
                component,
                gene,
                "yes",
                len(target_idx),
                safe_fmt(target_mean),
                safe_fmt(other_mean),
                safe_fmt(delta),
                support,
            ])
    return rows, component_gene_stats


def spatial_anchor():
    cols, raw_data = load_gse277104_qnorm(A01 / "GSE277104_qnorm_AOI_expression.tsv")
    a01_meta_rows = read_tsv(A01 / "GSE277104_AOI_metadata.tsv")
    series_meta = read_series_matrix_meta(A01 / "GSE277104_series_matrix.txt.gz")
    supplementary = series_meta.get("!Sample_supplementary_file_1", [[]])[0]
    dcc_to_meta = {}
    for dcc, meta in zip([dcc_name_from_supplementary(x) for x in supplementary], a01_meta_rows):
        dcc_to_meta[dcc] = meta
    matched_idx = [i for i, col in enumerate(cols) if col in dcc_to_meta]
    if not matched_idx:
        raise RuntimeError("No GSE277104 qnorm DCC columns matched GEO supplementary DCC names")
    cols = [cols[i] for i in matched_idx]
    data = {gene: [vals[i] for i in matched_idx] for gene, vals in raw_data.items()}
    meta_rows = [dcc_to_meta[col] for col in cols]
    aoi_types = [r["aoi_type"] for r in meta_rows]
    tumor_idx = [i for i, x in enumerate(aoi_types) if x == "Tumor"]
    vessel_idx = [i for i, x in enumerate(aoi_types) if x == "Vessel"]
    other_idx = [i for i, x in enumerate(aoi_types) if x not in {"Tumor", "Vessel"}]

    gene_rows = []
    module_rows = []
    spatial_gene_support = {}
    for component, genes in PREDEFINED_PROGRAMS.items():
        measured = [g for g in genes if g in data]
        if component == "malignant_hypoxia_adaptive":
            expected = "Tumor"
        elif component == "endothelial_angiogenesis_interaction":
            expected = "Vessel"
        else:
            expected = "Tumor_or_Vessel"
        module_tumor_scores = []
        module_vessel_scores = []
        for gene in measured:
            vals = data[gene]
            all_mean = mean(vals)
            all_sd = (sum((v - all_mean) ** 2 for v in vals) / max(len(vals) - 1, 1)) ** 0.5
            z = [(v - all_mean) / all_sd if all_sd else 0.0 for v in vals]
            tumor_mean = mean(vals[i] for i in tumor_idx)
            vessel_mean = mean(vals[i] for i in vessel_idx)
            other_mean = mean(vals[i] for i in other_idx)
            tumor_z = mean(z[i] for i in tumor_idx)
            vessel_z = mean(z[i] for i in vessel_idx)
            module_tumor_scores.append(tumor_z)
            module_vessel_scores.append(vessel_z)
            if expected == "Tumor":
                support = "spatial_expected_context" if tumor_z >= vessel_z else "spatial_measured_opposite_or_mixed"
            elif expected == "Vessel":
                support = "spatial_expected_context" if vessel_z >= tumor_z else "spatial_measured_opposite_or_mixed"
            else:
                support = "spatial_measured_context" if max(tumor_z, vessel_z) > -0.1 else "spatial_low_or_mixed"
            spatial_gene_support[(component, gene)] = support
            gene_rows.append([
                component,
                gene,
                "yes",
                expected,
                safe_fmt(tumor_mean),
                safe_fmt(vessel_mean),
                safe_fmt(other_mean),
                safe_fmt(tumor_z),
                safe_fmt(vessel_z),
                support,
            ])
        for gene in genes:
            if gene not in data:
                spatial_gene_support[(component, gene)] = "not_measured_in_GSE277104"
                gene_rows.append([component, gene, "no", expected, "", "", "", "", "", "not_measured_in_GSE277104"])
        module_rows.append([
            component,
            ";".join(measured),
            len(measured),
            expected,
            safe_fmt(mean(module_tumor_scores)),
            safe_fmt(mean(module_vessel_scores)),
            "spatial_anchor_available" if measured else "no_spatial_anchor",
        ])
    return gene_rows, module_rows, spatial_gene_support


def gene_universe(path, delimiter="\t"):
    out = set()
    with Path(path).open(encoding="utf-8-sig") as fh:
        next(fh)
        for line in fh:
            if line.strip():
                out.add(line.split(delimiter, 1)[0].strip())
    return out


def bulk_availability(final_genes):
    universes = {
        "GSE104580": gene_universe(A01 / "GSE104580_gene_expression.tsv"),
        "GSE14520_GPL3921": gene_universe(A01 / "GSE14520_GPL3921_gene_expression.tsv"),
        "GSE14520_GPL571": gene_universe(A01 / "GSE14520_GPL571_gene_expression.tsv"),
        "GSE277104": gene_universe(A01 / "GSE277104_qnorm_AOI_expression.tsv"),
    }
    rows = []
    for r in final_genes:
        gene = r["gene_symbol"]
        rows.append([
            gene,
            r["component"],
            *["yes" if gene in universes[k] else "no" for k in universes],
            "yes" if all(gene in universes[k] for k in universes) else "no",
        ])
    return rows


def main():
    A02.mkdir(parents=True, exist_ok=True)
    sc_rows_149, sc_stats_149 = single_cell_anchor(
        "GSE149614", "processed_GSE149614_candidate_gene_pseudobulk_counts.csv"
    )
    sc_rows_151, sc_stats_151 = single_cell_anchor(
        "GSE151530", "processed_GSE151530_HCC_candidate_gene_pseudobulk_counts.csv"
    )
    write_csv(
        "A02_single_cell_anchor_programs.csv",
        [
            "dataset",
            "component",
            "gene_symbol",
            "detected",
            "n_target_groups",
            "target_mean_log1p_count_per_cell",
            "other_mean_log1p_count_per_cell",
            "target_minus_other",
            "single_cell_support",
        ],
        sc_rows_149 + sc_rows_151,
    )

    spatial_gene_rows, spatial_module_rows, spatial_support = spatial_anchor()
    write_csv(
        "A02_spatial_anchor_programs.csv",
        [
            "component",
            "gene_symbol",
            "measured_in_GSE277104",
            "expected_spatial_context",
            "tumor_AOI_mean",
            "vessel_AOI_mean",
            "other_AOI_mean",
            "tumor_AOI_mean_z",
            "vessel_AOI_mean_z",
            "spatial_support",
        ],
        spatial_gene_rows,
    )
    write_csv(
        "A02_spatial_anchor_module_summary.csv",
        [
            "component",
            "measured_genes",
            "n_measured_genes",
            "expected_spatial_context",
            "module_tumor_mean_z",
            "module_vessel_mean_z",
            "module_spatial_status",
        ],
        spatial_module_rows,
    )

    transfer_universes = {
        "GSE104580": gene_universe(A01 / "GSE104580_gene_expression.tsv"),
        "GSE14520_GPL3921": gene_universe(A01 / "GSE14520_GPL3921_gene_expression.tsv"),
        "GSE14520_GPL571": gene_universe(A01 / "GSE14520_GPL571_gene_expression.tsv"),
        "GSE277104": gene_universe(A01 / "GSE277104_qnorm_AOI_expression.tsv"),
    }
    excluded_rows = []
    final_genes = []
    for component, genes in PREDEFINED_PROGRAMS.items():
        for gene in genes:
            sc149 = sc_stats_149.get(("GSE149614", component, gene), {}).get("support", "not_detected")
            sc151 = sc_stats_151.get(("GSE151530", component, gene), {}).get("support", "not_detected")
            spatial = spatial_support.get((component, gene), "not_measured_in_GSE277104")
            sc_supported = sc149 == "target_enriched" or sc151 == "target_enriched"
            spatial_measured = spatial != "not_measured_in_GSE277104"
            bulk_transferable = all(gene in u for u in transfer_universes.values())
            if sc_supported and spatial_measured:
                row = {
                    "ecotype_name": "Hypoxia_adaptive_Tumor_Endothelial_Macrophage_Ecotype_A02",
                    "component": component,
                    "gene_symbol": gene,
                    "inclusion_status": "A03_core",
                    "single_cell_anchor": ";".join([x for x in [f"GSE149614:{sc149}", f"GSE151530:{sc151}"] if x]),
                    "spatial_anchor": spatial,
                    "inclusion_rule": "predefined_gene;target_cell_state_enriched_in_scRNA;measured_in_GSE277104;available_in_GSE104580_GSE14520_GPL3921_GPL571_GSE277104",
                }
                if bulk_transferable:
                    final_genes.append(row)
                else:
                    excluded_rows.append([
                        component,
                        gene,
                        "excluded_after_bulk_transfer_check",
                        row["single_cell_anchor"],
                        spatial,
                        ";".join(k for k, u in transfer_universes.items() if gene not in u),
                    ])

    # Keep unique gene-component rows, preserving cross-component genes such as VEGFA in their biological component.
    write_csv(
        "A02_final_ecotype_gene_set.csv",
        [
            "ecotype_name",
            "component",
            "gene_symbol",
            "inclusion_status",
            "single_cell_anchor",
            "spatial_anchor",
            "inclusion_rule",
        ],
        [[r[k] for k in [
            "ecotype_name",
            "component",
            "gene_symbol",
            "inclusion_status",
            "single_cell_anchor",
            "spatial_anchor",
            "inclusion_rule",
        ]] for r in final_genes],
    )
    write_csv(
        "A02_excluded_candidate_genes.csv",
        [
            "component",
            "gene_symbol",
            "exclusion_reason",
            "single_cell_anchor",
            "spatial_anchor",
            "missing_transfer_dataset",
        ],
        excluded_rows,
    )

    availability_rows = bulk_availability(final_genes)
    write_csv(
        "A02_bulk_transfer_gene_availability.csv",
        [
            "gene_symbol",
            "component",
            "GSE104580",
            "GSE14520_GPL3921",
            "GSE14520_GPL571",
            "GSE277104",
            "available_in_all_transfer_datasets",
        ],
        availability_rows,
    )

    component_counts = defaultdict(int)
    for r in final_genes:
        component_counts[r["component"]] += 1
    all_available = sum(1 for r in availability_rows if r[-1] == "yes")
    verdict = "READY_FOR_A03" if final_genes and all_available == len(availability_rows) else "BLOCKED_NEED_REPAIR"

    summary = f"""# A02 ecotype construction summary

本阶段构建 single-cell / spatial anchored hypoxia-adaptive tumor-endothelial-macrophage ecotype。A02 未使用 GSE104580 response 标签筛选基因，未使用 GSE14520 survival/RFS/recurrence/treatment outcome 标签定义基因，未做差异分析、LASSO、随机森林、SVM、预测建模或正式生存分析。

## Biological rationale

本 ecotype 不是普通 TACE response signature，而是一个预治疗不适配状态的可转译框架：恶性肿瘤细胞缺氧适应、血管/内皮相互作用、SPP1/MMP9-like TAM 和基质重塑支持在空间上共同形成 tumor-vessel 相关生态位。后续 A03 只检验该预定义 ecotype 是否与 GSE104580 TACE non-response 相关。

## Single-cell anchor

- GSE149614 使用 candidate-gene pseudobulk，按 `Hepatocyte`、`Endothelial`、`Myeloid` 作为 tumor/endothelial/macrophage anchor。
- GSE151530 使用 HCC-filtered candidate-gene pseudobulk，按 `Malignant cells`、`TECs`、`TAMs` 作为 tumor/endothelial/macrophage anchor。
- 输出表：`A02_single_cell_anchor_programs.csv`。

## Spatial anchor

- GSE277104 使用 Tumor AOI 与 Vessel AOI，仅作为 biological spatial anchor。
- 输出 gene-level spatial anchor：`A02_spatial_anchor_programs.csv`。
- 输出 module-level spatial summary：`A02_spatial_anchor_module_summary.csv`。
- 本阶段不把 spatial anchor 写成 TACE 治疗证据，也不转向 MVI/perivascular 主线。

## Ecotype definition

最终 ecotype 名称：`Hypoxia_adaptive_Tumor_Endothelial_Macrophage_Ecotype_A02`。

纳入规则：

1. 基因必须来自 A01 预先指定 biological candidate programs。
2. 基因必须在 GSE149614 或 GSE151530 的目标细胞状态中呈 target-enriched。
3. 基因必须在 GSE277104 中可测，作为 spatial anchor 的可转译约束。
4. 基因必须通过 bulk-transfer availability 检查。
5. 不允许使用 GSE104580 response 标签或 GSE14520 outcome 标签增删基因。

最终 A03 core genes 数：{len(final_genes)}。

按 component 计数：

{chr(10).join(f'- {k}: {v}' for k, v in sorted(component_counts.items()))}

## 与普通 signature 的区别

- 不从 TACE response DEGs 出发。
- 不训练预测模型。
- 不把单一 hypoxia、angiogenesis 或 TAM signature 作为结论。
- 强制要求 tumor、endothelial/vessel、macrophage support 三个生态位 component 共同构成可转译模块。

## Bulk-transfer readiness

Gene availability 输出：`A02_bulk_transfer_gene_availability.csv`。

Final verdict: `{verdict}`。
"""
    (A02 / "A02_ecotype_construction_summary.md").write_text(summary, encoding="utf-8")

    decisions = """# A02 method decisions and forbidden actions

## Method decisions

- 使用 A01 预先指定候选 gene programs，不从 GSE104580 response 标签筛选基因。
- single-cell anchor 使用 GSE149614 与 GSE151530 HCC pseudobulk 的 target cell-state recovery。
- spatial anchor 使用 GSE277104 Tumor/Vessel AOI 的 gene/module expression context。
- A03 core ecotype genes 必须同时具备 single-cell target-state 支持、GSE277104 可测性和 bulk transfer 可测性。
- GSE104580/GSE14520 在 A02 只用于 gene availability，不参与 ecotype gene selection。

## Forbidden actions

- 不允许做 GSE104580 response DEGs 先行筛基因。
- 不允许 LASSO、随机森林、SVM 或其他预测建模。
- 不允许正式生存分析。
- 不允许把 OS/RFS 反客为主。
- 不允许把 GSE14520 写成 response validation cohort。
- 不允许根据显著性后验替换 endpoint。
- 不允许转向 ICI/TLS/MVI/perivascular/immune-exclusion/早复发主线。
- 不允许把本研究写成普通 TACE response gene signature。
"""
    (A02 / "A02_method_decisions_and_forbidden_actions.md").write_text(decisions, encoding="utf-8")

    final = f"""# final_A02_verdict

`{verdict}`

## 最终 ecotype / module / score 名称

`Hypoxia_adaptive_Tumor_Endothelial_Macrophage_Ecotype_A02`

## A03 可直接使用的 gene set 文件路径

`results/a02_ecotype_construction/A02_final_ecotype_gene_set.csv`

## A03 中允许做什么

- 将 A02 预定义 ecotype genes 转译到 GSE104580 bulk pretreatment expression。
- 计算 rank-based module/ecotype score、ssGSEA/AUCell/singscore 类 score。
- 检验该预定义 score 与 GSE104580 TACE non-response 的方向和强度。
- 仅作为后续链条的一步，保持 GSE14520 为 A04 TACE-treated outcome validation 和 resection-only specificity check。

## A03 中仍然禁止做什么

- 不允许做 GSE104580 response DEGs 后再回填 ecotype。
- 不允许 LASSO、随机森林、SVM 或其他预测建模。
- 不允许正式生存分析或把 OS/RFS 反客为主。
- 不允许把 GSE14520 写成 response validation cohort。
- 不允许根据显著性后验替换 endpoint 或增删 A02 gene set。
- 不允许转向 ICI/TLS/MVI/perivascular/immune-exclusion/早复发主线。
"""
    (A02 / "final_A02_verdict.md").write_text(final, encoding="utf-8")


if __name__ == "__main__":
    main()
