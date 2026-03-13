# TCGA-BRCA Breast Cancer Target Identification Pipeline
### A Multi-Phase Bioinformatics Workflow for Drug Target Discovery
  
> Integrating multi-source transcriptomic data to identify, validate, and prioritize druggable targets in breast cancer through differential expression, functional enrichment, TCGA validation, subtype analysis, survival analysis, and cell line confirmation.

---
## Project Overview

This pipeline integrates differential expression results from two independent breast cancer transcriptomic studies to identify a high-confidence set of tumor-upregulated genes. These genes are subsequently validated in the full TCGA-BRCA cohort, characterized by PAM50 molecular subtype, assessed for survival significance, confirmed in breast cancer cell lines, and scored for druggability — culminating in the selection of three priority drug targets for molecular docking.

### Core Question
*Which genes are consistently upregulated across independent breast cancer datasets, validated in large patient cohorts, subtype-specific, survival-relevant, expressed in cell lines, and structurally amenable to drug targeting?*

## Pipeline Summary

```
Project 1 DEGs (4,103 genes)  ──┐
                                 ├──► 666 Overlap Genes
Project 2 DEGs (1,328 genes)  ──┘         │
                                           ▼
                                  Functional Enrichment
                                  (GO, KEGG, Reactome)
                                           │
                                           ▼
                                  TCGA Validation (n=1,093 tumors)
                                  502 Confirmed Upregulated (93%)
                                           │
                                           ▼
                                  PAM50 Subtype Analysis
                                  (521 annotated samples)
                                           │
                                           ▼
                                  Survival Analysis
                                  (Cox regression, 502 genes)
                                           │
                                           ▼
                                  CCLE Cell Line Expression
                                  (depmap 22Q2)
                                           │
                                           ▼
                                  Evidence Scoring + Druggability Filter
                                           │
                                           ▼
                               Probable drug targets
                              
```

---

## Data Sources

| Source | Description | Access |
|--------|-------------|--------|
| **Project 1** | TCGA-BRCA 30 tumor vs 30 normal, TCGAbiolinks/DESeq2 | `TCGA_BRCA_HighConfidence_DEGs1_Project_1.csv` |
| **Project 2** | Salmon quantification of raw FASTQ data, DESeq2 | `WGS_Breast_cancer_Full_Differential_Expression_project_2.csv` |
| **TCGA-BRCA Full Cohort** | 1,093 tumor + 112 normal samples | `curatedTCGAData` Bioconductor (v2.0.1) |
| **CCLE/depmap** | Breast cancer cell line expression | `depmap` Bioconductor (22Q2) |

---

## Requirements

### R Version
R 4.5.2, Bioconductor 3.22

### Required Packages
```r
# CRAN
install.packages(c("dplyr", "tidyr", "ggplot2", "ggpubr", "pheatmap", 
                   "gridExtra", "tibble"))

# Bioconductor
BiocManager::install(c(
  "TCGAbiolinks", "DESeq2", "curatedTCGAData", "MultiAssayExperiment",
  "TCGAutils", "clusterProfiler", "enrichplot", "ReactomePA",
  "org.Hs.eg.db", "AnnotationDbi", "ggVennDiagram",
  "survival", "survminer", "depmap", "ExperimentHub"
))
```

### Important Notes
- Use `dplyr::select()` explicitly — masked by `AnnotationDbi`
- Use `stats::median()` explicitly — masked in depmap environment
- `curatedTCGAData` streams to local cache — requires internet on first run


## Phase 1 — Data Integrity & Cleaning

**Objective:** Load, standardize, and filter DEG results from both projects to produce clean upregulated gene lists.

**Key Steps:**
- Stripped Ensembl version suffixes (`ENSG00000xxx.21` → `ENSG00000xxx`)
- Filtered: `padj < 0.05` AND `log2FC > 1`
- Removed NAs and duplicates (kept lowest padj per symbol)

**Results:**
| Dataset | Total Genes | Upregulated (after filtering) |
|---------|-------------|-------------------------------|
| Project 1 | 9,383 | 4,103 |
| Project 2 | 29,161 | 1,328 |

---

## Phase 2 — Overlap Analysis

**Objective:** Identify genes consistently upregulated across both independent datasets.

**Method:** Gene symbol intersection using `intersect()` after deduplication.

**Result:** **666 overlap genes** (intersection of 4,103 and 1,328 upregulated gene sets)

---

<img width="2400" height="1800" alt="Venn_Diagram_Overlap" src="https://github.com/user-attachments/assets/59238ee6-a4c3-48f6-9031-cb0760610f18" />

> *Two-circle Venn diagram showing overlap between Project 1 (4,103 genes) and Project 2 (1,328 genes). The intersection of 666 genes is highlighted. This represents your high-confidence tumor-upregulated gene set carried forward for all downstream analyses.*

---

## Phase 3 — Functional Enrichment

**Objective:** Characterize the biological pathways and processes enriched among the 666 overlap genes.

**Methods:** GO (BP/MF/CC), KEGG, and Reactome enrichment via `clusterProfiler` and `ReactomePA`. Gene symbols converted to Entrez IDs via `bitr()` — 666/666 mapped successfully.

**Results:**
| Database | Significant Terms |
|----------|------------------|
| GO Biological Process | 390 |
| GO Molecular Function | 34 |
| GO Cellular Component | 64 |
| KEGG | 28 pathways |
| Reactome | 243 pathways |

**Dominant Biological Theme:** Cell Cycle Dysregulation & Chromosomal Instability
- Top GO-BP: nuclear division, chromosome segregation, mitotic sister chromatid segregation
- Top KEGG: Cell cycle, Fanconi anemia, Homologous recombination, p53 signaling
- Top Reactome: Cell Cycle Checkpoints, M Phase, DNA Replication, DNA methylation

<img width="1098" height="600" alt="GO_BP_Dotplot" src="https://github.com/user-attachments/assets/cc7053c4-1701-42e7-a987-1648867ea1e0" />

> *Dotplot of top 20 enriched GO Biological Process terms. Dot size = gene count, color = adjusted p-value.*

###  `Phase3_KEGG_Barplot.png`
> *Barplot of top KEGG pathways ranked by gene count. Cell cycle, Fanconi anemia, and homologous recombination pathways expected at top — consistent with chromosomal instability phenotype.*

<img width="893" height="488" alt="KEGG_Barplot" src="https://github.com/user-attachments/assets/4f560445-7e1d-4420-a951-9500584024ff" />

> *Dotplot of top Reactome pathways. Complements KEGG results with finer resolution — Cell Cycle Checkpoints, M Phase, and DNA Replication dominate.*

<img width="1098" height="600" alt="GO_BP_Emapplot" src="https://github.com/user-attachments/assets/bb7a9776-515c-49ba-9e70-117b524b7732" />
> *Enrichment map showing relationships between GO-BP terms as a network. Clustered nodes represent related biological processes*

## Phase 4 — TCGA Validation

**Objective:** Validate the 666 overlap genes in the full TCGA-BRCA cohort using independent expression data.

**Method:**
- Data streamed via `curatedTCGAData("BRCA", "RNASeq2GeneNorm", version="2.0.1")`
- Expression matrix: 20,501 genes × 1,212 samples
- Tumor samples (barcode pos 14–15 = "01"): 1,093
- Normal samples ("11"): 112
- Computed mean tumor/normal expression per gene
- log2FC = log2((mean_tumor + 1) / (mean_normal + 1))
- Confirmed upregulated: log2FC > 1

**Results:**
| Step | Count |
|------|-------|
| Genes submitted | 666 |
| Genes found in TCGA | 540 |
| Confirmed upregulated (log2FC > 1) | **502** |
| **Validation rate** | **93%** |

<img width="2400" height="2100" alt="Tumor_vs_Normal_Scatter" src="https://github.com/user-attachments/assets/b753ff1c-36b0-44b8-8f53-24d7e8259ab1" />

> *Scatter plot of mean tumor expression vs mean normal expression for all genes. Points above the diagonal represent upregulated genes. Validated genes (log2FC > 1) highlighted.*

## Phase 5 — Subtype Analysis

**Objective:** Determine which PAM50 molecular subtypes show highest expression of the 502 validated genes.

**Method:**
- PAM50 annotations extracted from `colData(brca_mae)` via `PAM50.mRNA` column
- 521 tumor samples matched to subtype annotation
- Mean expression computed per gene per subtype
- Top subtype identified per gene via `slice_max()`

**Subtype Sample Distribution:**
| Subtype | Samples |
|---------|---------|
| Luminal A | 231 |
| Luminal B | 127 |
| Basal-like | 97 |
| HER2-enriched | 58 |
| Normal-like | 8 |

**Gene Enrichment per Subtype (top subtype per gene):**
| Subtype | Genes | % of 502 |
|---------|-------|----------|
| Basal-like | 221 | 44% |
| HER2-enriched | 113 | 22% |
| Luminal B | 63 | 13% |
| Normal-like | 56 | 11% |
| Luminal A | 49 | 10% |

> Despite Luminal A having the most patients (231), Basal-like dominates gene enrichment (221 genes) — reflecting the transcriptionally extreme nature of triple-negative/Basal-like tumors. The 502-gene signature is enriched in the two most aggressive subtypes (Basal + HER2 = 66%).

<img width="2400" height="1800" alt="Subtype_Distribution" src="https://github.com/user-attachments/assets/16e7082b-58c6-4fce-9cf3-bdcef8b0fc46" />
> *Bar chart showing number of validated genes with highest expression per PAM50 subtype. Basal-like bar expected to be tallest (221 genes), followed by HER2-enriched (113). Color-coded by subtype using standard PAM50 colors (red=Basal, orange=HER2, green=LumA, blue=LumB, purple=Normal-like).*

<img width="4200" height="3600" alt="Top20_Subtype_Boxplots" src="https://github.com/user-attachments/assets/ce87a3ad-f521-43b3-9ac9-c705ac996d93" />

> *Faceted boxplot grid (4×5) showing expression of top 20 validated genes across all 5 PAM50 subtypes. Each panel = one gene, x-axis = subtype, y-axis = normalized expression. Three patterns visible: (1) Basal-like enriched genes (STAT1, S100A11), (2) HER2-enriched dominant (ERBB2 — positive control), (3) Pan-subtype stromal genes (COL1A1, FN1, POSTN).*


## Phase 6 — Survival Analysis

**Objective:** Identify which of the 502 validated genes show statistically significant association with overall survival in TCGA-BRCA patients.

**Method:**
- Survival data: `vital_status` (0/1), `days_to_death`, `days_to_last_followup` from `colData(brca_mae)`
- OS time = days_to_death if dead, else days_to_last_followup
- KM curves: median expression split (High vs Low) per gene
- Cox regression: continuous expression, all 502 genes, BH multiple testing correction

**Cohort:**
| Metric | Value |
|--------|-------|
| Total samples | 1,079 |
| Death events | 152 |
| Censored (alive) | 927 |
| Event rate | 14% |

> Low event rate is a known characteristic of TCGA-BRCA — an early-stage cohort with ~85% 5-year survival. This limits statistical power for survival analysis.

**Cox Regression Results:**
- Genes tested: 502
- Significant after BH correction (padj < 0.05): **4 genes**
- All 4 had HR > 1 (worse survival with high expression)

| Gene | HR | padj | Interpretation |
|------|----|------|----------------|
| MAL2 | ~1.0001 | < 0.05 | Membrane trafficking, breast cancer progression |
| MMP13 | ~1.0008 | < 0.05 | ECM degradation, invasion |
| SHCBP1 | ~1.001 | < 0.05 | Spindle checkpoint, mitotic driver |
| ULBP2 | ~1.004 | < 0.05 | NK cell ligand, immune evasion |

## Phase 7 — Cell Line Expression (CCLE)

**Objective:** Confirm that candidate genes are expressed in breast cancer cell lines (not just patient tumor tissue), establishing experimental tractability for in vitro validation and drug testing.

**Method:**
- Data: `depmap` Bioconductor package, TPM_22Q2 (EH7556) + metadata_22Q2 (EH7558)
- Filtered to breast cancer cell lines via `primary_disease == "Breast Cancer"`
- Expression column: `rna_expression` (TPM)

**Results:**
- Breast cancer cell lines identified: confirmed panel
- Validated genes found in depmap: confirmed from 502

<img width="1000" height="700" alt="Top30_CCLE_Expression" src="https://github.com/user-attachments/assets/af9112cd-90fb-4d7d-a691-22aaa175eab8" />

> *Horizontal bar chart of top 30 validated genes by mean TPM in breast cancer cell lines. Color gradient from light orange (low) to dark red (high). Shows which of the 502 validated genes are most robustly expressed in experimental cell line models.*

## Phase 8 — Target Prioritization Scoring

**Objective:** Objectively rank all 502 validated genes using a multi-criteria evidence scoring system to select the top probable drug targets.

### Scoring Criteria (Maximum 14 points)

| Criterion | Max Points | Rationale |
|-----------|-----------|-----------|
| Mean log2FC (P1+P2) | 3 | Expression strength across discovery datasets |
| TCGA log2FC | 3 | Validation strength in large cohort |
| Subtype specificity | 2 | Aggressive subtype (Basal/HER2) enrichment |
| Survival significance | 3 | Cox regression padj and HR direction |
| CCLE expression | 3 | Cell line tractability |
| Druggability bonus | +2 | Known oncogenic/druggable gene families |
| Stromal penalty | -3 | ECM/stromal genes deprioritized |


<img width="1500" height="1200" alt="Phase8_Target_Prioritization_Scores" src="https://github.com/user-attachments/assets/d8c3ff3e-b432-4090-be29-e73ca9e090a7" />

> *Horizontal bar chart of top 20 scored genes. Bars color-coded by top PAM50 subtype (red=Basal-like, orange=HER2-enriched, blue=Luminal B). Score labels shown at end of each bar. PLK1, TOP2A, CDK1 visible in upper portion of chart.*

## Citation & Acknowledgements

Data sources:
- TCGA Research Network: https://www.cancer.gov/tcga
- curatedTCGAData: Ramos et al., Cancer Research (2020)
- DepMap: Broad Institute Cancer Dependency Map
- Reactome: Fabregat et al., Nucleic Acids Research (2018)
- KEGG: Kanehisa & Goto, Nucleic Acids Research (2000)

*Pipeline developed in R 4.5.2 | Bioconductor 3.22 |*
