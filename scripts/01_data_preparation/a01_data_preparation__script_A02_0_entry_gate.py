import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path


BASE = Path(__file__).resolve().parent
LOG = BASE / "log_A02_0_entry_gate.txt"


def log(msg):
    print(msg)
    with LOG.open("a", encoding="utf-8") as fh:
        fh.write(str(msg) + "\n")


def open_text(path):
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def write_csv(path, header, rows):
    with (BASE / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def write_tsv(path, header, rows):
    with (BASE / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(header)
        w.writerows(rows)


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
            parts = [x.strip('"') for x in line.split("\t")]
            if in_table:
                if header is None:
                    header = parts
                else:
                    rows.append(parts)
            elif line.startswith("!"):
                meta[parts[0]].append(parts[1:])
    return meta, header, rows


def load_probe_map(path):
    out = {}
    with (BASE / path).open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            probe = r.get("probe_id", "")
            gene = r.get("gene_symbol", "")
            if probe and gene:
                out[probe] = gene
    return out


def load_candidate_genes():
    genes = []
    seen = set()
    with (BASE / "candidate_gene_availability.tsv").open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            g = r["gene"]
            if g not in seen:
                genes.append(g)
                seen.add(g)
    return genes


def platform_audit():
    probe_map = load_probe_map("GPL570_probe_to_gene.tsv")
    clinical = {}
    with (BASE / "GSE14520_Fako_mapped_clinical_metadata.tsv").open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if r["gsm"]:
                clinical[r["gsm"]] = r

    audit_rows = []
    sample_to_platform = {}
    for platform, matrix_file, gene_file in [
        ("GPL3921", "GSE14520-GPL3921_series_matrix.txt.gz", "GSE14520_GPL3921_gene_expression.tsv"),
        ("GPL571", "GSE14520-GPL571_series_matrix.txt.gz", "GSE14520_GPL571_gene_expression.tsv"),
    ]:
        meta, header, rows = read_series_matrix(BASE / matrix_file)
        samples = [x.strip('"') for x in header[1:]]
        for s in samples:
            sample_to_platform[s] = platform
        probes = [r[0].strip('"') for r in rows]
        mapped_probes = [p for p in probes if p in probe_map]
        genes = set()
        with (BASE / gene_file).open(encoding="utf-8") as fh:
            next(fh)
            for line in fh:
                if line.strip():
                    genes.add(line.split("\t", 1)[0])
        clinical_on_platform = [clinical[s] for s in samples if s in clinical]
        tace = sum(1 for r in clinical_on_platform if r["analysis_context"] == "TACE_treated")
        resection = sum(1 for r in clinical_on_platform if r["analysis_context"] == "resection_only")
        other = sum(1 for r in clinical_on_platform if r["analysis_context"] == "other_therapy")
        missing = sum(1 for r in clinical_on_platform if r["analysis_context"] == "missing_survival_data")
        audit_rows.append([
            platform,
            matrix_file,
            len(samples),
            len(probes),
            len(mapped_probes),
            f"{len(mapped_probes) / max(len(probes), 1):.4f}",
            len(genes),
            sum(1 for s in samples if s in clinical),
            tace,
            resection,
            other,
            missing,
            "GPL570_probe_table_conservative_mapping",
            "usable_for_module_score_transfer_with_caveat",
            "Use only for TACE-treated post-TACE outcome validation and resection-only specificity; not response validation.",
        ])

    therapy_counts = Counter()
    platform_counts = Counter()
    for gsm, r in clinical.items():
        therapy_counts[r["analysis_context"]] += 1
        platform_counts[(r["analysis_context"], sample_to_platform.get(gsm, "not_in_matrix"))] += 1

    write_csv(
        "GSE14520_platform_mapping_audit.csv",
        [
            "platform",
            "matrix_file",
            "n_matrix_samples",
            "n_matrix_probes",
            "n_probes_mapped_by_current_probe_table",
            "probe_mapping_fraction",
            "n_gene_symbols_in_gene_matrix",
            "n_fako_clinical_samples_on_platform",
            "n_TACE_treated_on_platform",
            "n_resection_only_on_platform",
            "n_other_therapy_on_platform",
            "n_missing_survival_on_platform",
            "mapping_source_used",
            "entry_gate_conclusion",
            "allowed_downstream_use",
        ],
        audit_rows,
    )
    return audit_rows, therapy_counts, platform_counts


def prepare_gse149614(candidate_genes):
    log("Preparing GSE149614 candidate-gene pseudobulk")
    meta = {}
    sample_set = set()
    patient_set = set()
    celltype_counts = Counter()
    with (BASE / "GSE149614_cell_metadata.tsv").open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fields = reader.fieldnames or []
        for r in reader:
            cell = r["Cell"]
            group = f'{r["sample"]}|{r["site"]}|{r["celltype"]}'
            meta[cell] = {
                "group": group,
                "sample": r["sample"],
                "site": r["site"],
                "patient": r["patient"],
                "celltype": r["celltype"],
            }
            sample_set.add(r["sample"])
            patient_set.add(r["patient"])
            celltype_counts[r["celltype"]] += 1

    group_cells = Counter(v["group"] for v in meta.values())
    group_keys = sorted(group_cells)
    group_index = {g: i for i, g in enumerate(group_keys)}
    write_tsv(
        "metadata_GSE149614_pseudobulk_group_metadata.tsv",
        ["group_id", "sample", "site", "celltype", "n_cells"],
        [[g, *g.split("|"), group_cells[g]] for g in group_keys],
    )

    candidate_set = set(candidate_genes)
    sums = {g: [0.0] * len(group_keys) for g in candidate_genes}
    detected = set()
    n_genes_raw = 0
    with open_text(BASE / "GSE149614_HCC.scRNAseq.S71915.count.txt.gz") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        col_group = [group_index[meta[c]["group"]] if c in meta else None for c in header]
        for line in fh:
            if not line.strip():
                continue
            n_genes_raw += 1
            gene, values = line.rstrip("\n").split("\t", 1)
            if gene not in candidate_set:
                continue
            detected.add(gene)
            vals = values.split("\t")
            arr = sums[gene]
            for i, x in enumerate(vals):
                gi = col_group[i]
                if gi is not None and x != "0":
                    arr[gi] += float(x)

    write_csv(
        "processed_GSE149614_candidate_gene_pseudobulk_counts.csv",
        ["gene_symbol"] + group_keys,
        [[g] + sums[g] for g in candidate_genes if g in detected],
    )
    return {
        "dataset": "GSE149614",
        "object_type": "candidate_gene_pseudobulk_counts",
        "object_path": str(BASE / "processed_GSE149614_candidate_gene_pseudobulk_counts.csv"),
        "group_metadata_path": str(BASE / "metadata_GSE149614_pseudobulk_group_metadata.tsv"),
        "n_cells": len(meta),
        "n_samples": len(sample_set),
        "n_patients": len(patient_set),
        "n_genes_raw": n_genes_raw,
        "n_candidate_genes_detected": len(detected),
        "hcc_filter_status": "GSE149614 HCC dataset; retained Tumor/PVTT/Normal/Lymph sites with site metadata",
        "celltype_field": "celltype",
        "major_celltypes": ";".join(f"{k}:{v}" for k, v in sorted(celltype_counts.items())),
        "qc_status": "ready_for_A02_cell_state_recovery_input",
    }


def prepare_gse151530(candidate_genes):
    log("Preparing GSE151530 HCC-filtered candidate-gene pseudobulk")
    hcc_samples = set()
    with (BASE / "GSE151530_sample_cancer_type.tsv").open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if "Hepatocellular" in r["cancer_type"]:
                hcc_samples.add(r["S_ID"])

    meta_by_cell = {}
    sample_set = set()
    celltype_counts = Counter()
    with (BASE / "GSE151530_cell_metadata.tsv").open(encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for r in reader:
            if r["S_ID"] not in hcc_samples:
                continue
            group = f'{r["S_ID"]}|{r["Type"]}'
            meta_by_cell[r["Cell"]] = {"sample": r["S_ID"], "type": r["Type"], "group": group}
            sample_set.add(r["S_ID"])
            celltype_counts[r["Type"]] += 1

    group_cells = Counter(v["group"] for v in meta_by_cell.values())
    group_keys = sorted(group_cells)
    group_index = {g: i for i, g in enumerate(group_keys)}
    write_tsv(
        "metadata_GSE151530_HCC_pseudobulk_group_metadata.tsv",
        ["group_id", "S_ID", "celltype", "n_cells"],
        [[g, *g.split("|"), group_cells[g]] for g in group_keys],
    )

    genes = []
    with open_text(BASE / "GSE151530_genes.tsv.gz") as fh:
        for line in fh:
            p = line.rstrip("\n").split("\t")
            genes.append(p[-1])
    candidate_set = set(candidate_genes)
    gene_index_to_symbol = {i + 1: g for i, g in enumerate(genes) if g in candidate_set}

    barcodes = []
    with open_text(BASE / "GSE151530_barcodes.tsv.gz") as fh:
        for line in fh:
            barcodes.append(line.strip())
    col_to_group = {}
    for i, bc in enumerate(barcodes, start=1):
        if bc in meta_by_cell:
            col_to_group[i] = group_index[meta_by_cell[bc]["group"]]

    sums = {g: [0.0] * len(group_keys) for g in candidate_genes}
    detected = set()
    mtx_shape = ""
    with open_text(BASE / "GSE151530_matrix.mtx.gz") as fh:
        for line in fh:
            if line.startswith("%"):
                continue
            mtx_shape = line.strip()
            break
        for line in fh:
            r, c, v = line.split()
            gene = gene_index_to_symbol.get(int(r))
            if not gene:
                continue
            gi = col_to_group.get(int(c))
            if gi is None:
                continue
            sums[gene][gi] += float(v)
            detected.add(gene)

    write_csv(
        "processed_GSE151530_HCC_candidate_gene_pseudobulk_counts.csv",
        ["gene_symbol"] + group_keys,
        [[g] + sums[g] for g in candidate_genes if g in detected],
    )
    return {
        "dataset": "GSE151530",
        "object_type": "HCC_filtered_candidate_gene_pseudobulk_counts",
        "object_path": str(BASE / "processed_GSE151530_HCC_candidate_gene_pseudobulk_counts.csv"),
        "group_metadata_path": str(BASE / "metadata_GSE151530_HCC_pseudobulk_group_metadata.tsv"),
        "n_cells": len(meta_by_cell),
        "n_samples": len(sample_set),
        "n_patients": "",
        "n_genes_raw": len(genes),
        "n_candidate_genes_detected": len(detected),
        "hcc_filter_status": f"HCC-filtered by sample_cancer_type; HCC_samples={len(hcc_samples)}; matrix_shape={mtx_shape}",
        "celltype_field": "Type",
        "major_celltypes": ";".join(f"{k}:{v}" for k, v in sorted(celltype_counts.items())),
        "qc_status": "ready_for_A02_cell_state_recovery_input",
    }


def main():
    LOG.write_text("", encoding="utf-8")
    candidate_genes = load_candidate_genes()
    audit_rows, therapy_counts, platform_counts = platform_audit()
    sc_rows = [prepare_gse149614(candidate_genes), prepare_gse151530(candidate_genes)]
    write_csv(
        "scRNA_object_preparation_summary.csv",
        [
            "dataset",
            "object_type",
            "object_path",
            "group_metadata_path",
            "n_cells",
            "n_samples",
            "n_patients",
            "n_genes_raw",
            "n_candidate_genes_detected",
            "hcc_filter_status",
            "celltype_field",
            "major_celltypes",
            "qc_status",
        ],
        [[r[k] for k in [
            "dataset",
            "object_type",
            "object_path",
            "group_metadata_path",
            "n_cells",
            "n_samples",
            "n_patients",
            "n_genes_raw",
            "n_candidate_genes_detected",
            "hcc_filter_status",
            "celltype_field",
            "major_celltypes",
            "qc_status",
        ]] for r in sc_rows],
    )

    gse14520_ready = all(float(r[5]) > 0.95 and int(r[8]) + int(r[9]) > 0 for r in audit_rows)
    sc_ready = all(r["n_candidate_genes_detected"] >= 30 and r["n_cells"] > 0 for r in sc_rows)
    verdict = "READY_FOR_A02" if gse14520_ready and sc_ready else "BLOCKED_NEED_REPAIR"

    summary = f"""# A02.0 entry gate summary

本轮只补齐 A01 verdict 中进入 A02 前的两个条件：GSE14520 平台映射可用性，以及 GSE149614/GSE151530 单细胞对象化输入。未进行 GSE104580 response 差异分析、LASSO、预测建模或正式生存分析。

## 1. GSE14520 / Fako 平台映射审计

- 当前 GSE14520 表达矩阵包含 GPL3921 与 GPL571 两个 Series Matrix。
- A01 平铺流程中使用 GPL570 probe-to-gene 表进行保守映射；A02.0 对两个矩阵 probe ID 的可映射比例进行了重新审计。
- 平台审计表：`GSE14520_platform_mapping_audit.csv`。
- 映射结论：当前映射足以支持后续 rank/module/ecotype score 的转译，但需要在方法中写明 probe annotation 使用保守 Affymetrix probe symbol 映射，不能把平台整合写成无批次问题。
- 允许用途：A04 中作为 TACE-treated post-TACE outcome validation 与 resection-only specificity check。
- 禁止用途：不得写成 TACE response validation cohort。

## 2. Single-cell 对象准备

- GSE149614 已生成候选基因 pseudobulk count 输入：`processed_GSE149614_candidate_gene_pseudobulk_counts.csv`。
- GSE149614 group metadata：`metadata_GSE149614_pseudobulk_group_metadata.tsv`。
- GSE151530 已按 HCC 样本过滤，生成候选基因 pseudobulk count 输入：`processed_GSE151530_HCC_candidate_gene_pseudobulk_counts.csv`。
- GSE151530 group metadata：`metadata_GSE151530_HCC_pseudobulk_group_metadata.tsv`。
- 对象级 QC 汇总：`scRNA_object_preparation_summary.csv`。
- 候选基因来自 A01 预先指定集合，不来自 GSE104580 response 标签。

## Entry Gate 结论

`{verdict}`

如果进入 A02，正式分析只能使用本轮列出的对象路径，并继续禁止差异分析先行、LASSO、预测模型和正式生存分析。
"""
    (BASE / "A02_0_entry_gate_summary.md").write_text(summary, encoding="utf-8")

    allowed = "\n".join([
        "- GSE104580 bulk matrix: `GSE104580_gene_expression.tsv`",
        "- GSE104580 metadata: `GSE104580_sample_metadata.tsv`",
        "- GSE14520 GPL3921 matrix: `GSE14520_GPL3921_gene_expression.tsv`",
        "- GSE14520 GPL571 matrix: `GSE14520_GPL571_gene_expression.tsv`",
        "- GSE14520 Fako metadata: `GSE14520_Fako_mapped_clinical_metadata.tsv`",
        "- GSE277104 AOI matrix: `GSE277104_qnorm_AOI_expression.tsv`",
        "- GSE277104 AOI metadata: `GSE277104_AOI_metadata.tsv`",
        "- GSE149614 pseudobulk: `processed_GSE149614_candidate_gene_pseudobulk_counts.csv`",
        "- GSE149614 pseudobulk metadata: `metadata_GSE149614_pseudobulk_group_metadata.tsv`",
        "- GSE151530 HCC pseudobulk: `processed_GSE151530_HCC_candidate_gene_pseudobulk_counts.csv`",
        "- GSE151530 HCC pseudobulk metadata: `metadata_GSE151530_HCC_pseudobulk_group_metadata.tsv`",
    ])
    final = f"""# final_A02_0_verdict

`{verdict}`

## 允许用于 A02 的数据对象

{allowed}

## A02 禁止事项

- 不允许做 GSE104580 response DEGs 先行筛基因。
- 不允许 LASSO、随机森林、SVM 或其他预测建模。
- 不允许正式生存分析或把 OS/RFS 反客为主。
- 不允许把 GSE14520 写成 response validation cohort。
- 不允许根据显著性后验替换 endpoint。
- 不允许转向 ICI/TLS/MVI/perivascular/immune-exclusion/早复发主线。

## 限制

- GSE14520 平台映射可用于 module/ecotype score 转译，但需保留 Affymetrix probe annotation 保守映射说明。
- 单细胞对象为候选基因 pseudobulk/cell-state recovery 输入，不是正式 Seurat/AnnData 全量对象；A02 中如需全量对象，可在不使用 response 标签筛基因的前提下另行转换。
"""
    (BASE / "final_A02_0_verdict.md").write_text(final, encoding="utf-8")
    log(f"A02.0 completed with verdict {verdict}")


if __name__ == "__main__":
    main()
