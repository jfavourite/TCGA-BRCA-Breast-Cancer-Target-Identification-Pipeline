# TCGA-BRCA Breast Cancer Target Identification Pipeline
### A Multi-Phase Bioinformatics Workflow for Drug Target Discovery

> **PhD Research Project**  
> Integrating multi-source transcriptomic data to identify, validate, and prioritize druggable targets in breast cancer through differential expression, functional enrichment, TCGA validation, subtype analysis, survival analysis, and cell line confirmation.

---

## Table of Contents
1. [Project Overview](#project-overview)
2. [Pipeline Summary](#pipeline-summary)
3. [Data Sources](#data-sources)
4. [Requirements](#requirements)
5. [Phase 1 — Data Integrity & Cleaning](#phase-1--data-integrity--cleaning)
6. [Phase 2 — Overlap Analysis](#phase-2--overlap-analysis)
7. [Phase 3 — Functional Enrichment](#phase-3--functional-enrichment)
8. [Phase 4 — TCGA Validation](#phase-4--tcga-validation)
9. [Phase 5 — Subtype Analysis](#phase-5--subtype-analysis)
10. [Phase 6 — Survival Analysis](#phase-6--survival-analysis)
11. [Phase 7 — Cell Line Expression (CCLE)](#phase-7--cell-line-expression-ccle)
12. [Phase 8 — Target Prioritization Scoring](#phase-8--target-prioritization-scoring)
13. [Phase 9 — Druggability & Protein Structure](#phase-9--druggability--protein-structure)
14. [Final Selected Targets](#final-selected-targets)
15. [Output Files](#output-files)
16. [Key Results Summary](#key-results-summary)

---

## Project Overview

This pipeline integrates differential expression results from two independent breast cancer transcriptomic studies to identify a high-confidence set of tumor-upregulated genes. These genes are subsequently validated in the full TCGA-BRCA cohort, characterized by PAM50 molecular subtype, assessed for survival significance, confirmed in breast cancer cell lines, and scored for druggability — culminating in the selection of three priority drug targets for molecular docking.

### Core Question
*Which genes are consistently upregulated across independent breast cancer datasets, validated in large patient cohorts, subtype-specific, survival-relevant, expressed in cell lines, and structurally amenable to drug targeting?*

### Final Targets Selected
| Target | Subtype | Rationale |
|--------|---------|-----------|
| **PLK1** | Basal-like | Mitotic kinase, deep ATP pocket, clinical inhibitors available |
| **TOP2A** | HER2-enriched | Clinically validated, highest CCLE expression, FDA-approved inhibitors |
| **CDK1** | Basal-like | Master mitotic regulator, cannot be bypassed, excellent PDB structures |

---

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
                              PLK1 ── TOP2A ── CDK1
                              (Final Docking Candidates)
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
- Always use `dplyr::select()` explicitly — masked by `AnnotationDbi`
- Always use `stats::median()` explicitly — masked in depmap environment
- `curatedTCGAData` streams to local cache — requires internet on first run
- Cache location (Windows): `C:\Users\User\AppData\Local/R/cache/R/ExperimentHub`

---

## Phase 1 — Data Integrity & Cleaning

**Objective:** Load, standardize, and filter DEG results from both projects to produce clean upregulated gene lists.

**Key Steps:**
- Standardized column names (`Symbol` → `symbol`)
- Stripped Ensembl version suffixes (`ENSG00000xxx.21` → `ENSG00000xxx`)
- Filtered: `padj < 0.05` AND `log2FC > 1`
- Removed NAs and duplicates (kept lowest padj per symbol)

**Results:**
| Dataset | Total Genes | Upregulated (after filtering) |
|---------|-------------|-------------------------------|
| Project 1 | 9,383 | 4,103 |
| Project 2 | 29,161 | 1,328 |

**Output Files:**
- `Phase1_Project1_Upregulated_Clean.csv`
- `Phase1_Project2_Upregulated_Clean.csv`

---

## Phase 2 — Overlap Analysis

**Objective:** Identify genes consistently upregulated across both independent datasets.

**Method:** Gene symbol intersection using `intersect()` after deduplication.

**Result:** **666 overlap genes** (intersection of 4,103 and 1,328 upregulated gene sets)

---

###  `Phase2_Venn_Diagram_Overlap.png`
> *Two-circle Venn diagram showing overlap between Project 1 (4,103 genes) and Project 2 (1,328 genes). The intersection of 666 genes is highlighted. This represents your high-confidence tumor-upregulated gene set carried forward for all downstream analyses.*

---

**Output Files:**
- `Phase2_Venn_Diagram_Overlap.png`
- `Overlap_Genes.csv` — 666 genes with log2FC and padj from both projects

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

**Three Actionable Gene Clusters Identified:**
1. **Mitotic kinases** (CDK1, AURKB, PLK1, BUB1) — druggable, crystal structures available
2. **DNA repair genes** (BRCA2, RAD51, FANCD2) — synthetic lethality context
3. **Chromatin remodeling** — better as biomarkers

---

###  `Phase3_GO_BP_Dotplot.png`
> *Dotplot of top 20 enriched GO Biological Process terms. Dot size = gene count, color = adjusted p-value. Expected to show strong enrichment for mitotic/cell cycle terms confirming oncogenic transcriptional program.*

###  `Phase3_KEGG_Barplot.png`
> *Barplot of top KEGG pathways ranked by gene count. Cell cycle, Fanconi anemia, and homologous recombination pathways expected at top — consistent with chromosomal instability phenotype.*

### `Phase3_Reactome_Dotplot.png`
> *Dotplot of top Reactome pathways. Complements KEGG results with finer resolution — Cell Cycle Checkpoints, M Phase, and DNA Replication dominate.*

###  `Phase3_GO_BP_Emapplot.png`
> *Enrichment map showing relationships between GO-BP terms as a network. Clustered nodes represent related biological processes — expect to see a large cell cycle/mitosis cluster and a separate DNA repair cluster.*

---

**Output Files:**
- `Phase3_GO_BP_Results.csv`, `Phase3_GO_MF_Results.csv`, `Phase3_GO_CC_Results.csv`
- `Phase3_KEGG_Results.csv`, `Phase3_Reactome_Results.csv`
- `Phase3_GO_BP_Dotplot.png`, `Phase3_GO_MF_Dotplot.png`
- `Phase3_KEGG_Barplot.png`, `Phase3_Reactome_Dotplot.png`
- `Phase3_GO_BP_Emapplot.png`

---

## Phase 4 — TCGA Validation

**Objective:** Validate the 666 overlap genes in the full TCGA-BRCA cohort (700+ patients) using independent expression data.

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

> A validation rate of 93% is publication-worthy. Rates above 70% are considered strong independent validation.

---

###  `Phase4_TCGA_log2FC_Distribution.png`
> *Histogram of TCGA log2FC values for all 540 found genes. Should show a right-skewed distribution with the majority of genes having log2FC > 1, confirming broad upregulation in the full cohort.*

###  `Phase4_Tumor_vs_Normal_Scatter.png`
> *Scatter plot of mean tumor expression vs mean normal expression for all genes. Points above the diagonal represent upregulated genes. Validated genes (log2FC > 1) highlighted — expected to show clear separation from normal tissue.*

### `Phase4_Top30_Validated_Genes.png`
> *Bar chart of top 30 genes by TCGA log2FC. Genes like COL10A1, IBSP, MMP1 expected at top. Confirms the most dramatically overexpressed genes in 1,093 TCGA tumor samples.*

---

**Output Files:**
- `Phase4_All_Genes_TCGA_Validation.csv`
- `Phase4_Overlap_TCGA_Validated_Full.csv`
- `Phase4_Confirmed_Upregulated_Genes.csv` — 502 final validated genes

---

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

---

###  `Phase5_Subtype_Distribution.png`
> *Bar chart showing number of validated genes with highest expression per PAM50 subtype. Basal-like bar expected to be tallest (221 genes), followed by HER2-enriched (113). Color-coded by subtype using standard PAM50 colors (red=Basal, orange=HER2, green=LumA, blue=LumB, purple=Normal-like).*

###  `Phase5_Top20_Subtype_Boxplots.png`
> *Faceted boxplot grid (4×5) showing expression of top 20 validated genes across all 5 PAM50 subtypes. Each panel = one gene, x-axis = subtype, y-axis = normalized expression. Three patterns visible: (1) Basal-like enriched genes (STAT1, S100A11), (2) HER2-enriched dominant (ERBB2 — positive control), (3) Pan-subtype stromal genes (COL1A1, FN1, POSTN).*

###  `Phase5_Subtype_Heatmap.png` or `heatmap_top40_subtypes.png`
> *Hierarchical clustering heatmap of top 40–50 validated genes × 5 subtypes. Z-score scaled rows. Column clustering groups HER2-enriched with Luminal A, and Basal-like with Luminal B — reflecting known transcriptional relationships. Five gene clusters visible: (1) Normal-like enriched (SPP1, CRABP2), (2) HER2-dominant (ERBB2, CEACAM6), (3) Stromal/ECM pan-subtype (COL genes, FN1), (4) Basal-like specific (PARP1, KPNA2, STMN1, STAT1), (5) HER2+LumB shared (SDC1).*

---

**Output Files:**
- `Phase5_Subtype_Mean_Expression.csv`
- `Phase5_Gene_Top_Subtype.csv`
- `Phase5_Subtype_Enrichment_Count.csv`

---

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

---

### `Phase6_KM_PLK1.png`
> *Kaplan-Meier survival curve for PLK1 — High (red) vs Low (blue) expression groups. X-axis = days, Y-axis = survival probability. Tick marks = censored patients. Include log-rank p-value. Note: p may be non-significant due to low event numbers — report this as a cohort limitation.*

### `Phase6_KM_TOP2A.png`
> *KM curve for TOP2A expression vs overall survival. Same format as PLK1 curve.*

###  `Phase6_KM_CDK1.png`
> *KM curve for CDK1 expression vs overall survival.*

### `Phase6_Forest_Plot.png`
> *Forest plot of the 4 Cox-significant genes (MAL2, MMP13, SHCBP1, ULBP2). Each gene shown as a point estimate (HR) with 95% confidence interval. Dashed vertical line at HR=1. All points to the right of 1 indicating risk genes. Note the compressed x-axis scale (1.000–1.006) reflecting modest effect sizes in this low-event cohort.*

---

**Output Files:**
- `Phase6_Cox_Survival_Results.csv` — Cox results for all 502 genes
- `Phase6_KM_[GENE].png` — Individual KM curves for 8 candidate genes
- `Phase6_Forest_Plot.png`

---

## Phase 7 — Cell Line Expression (CCLE)

**Objective:** Confirm that candidate genes are expressed in breast cancer cell lines (not just patient tumor tissue), establishing experimental tractability for in vitro validation and drug testing.

**Method:**
- Data: `depmap` Bioconductor package, TPM_22Q2 (EH7556) + metadata_22Q2 (EH7558)
- Filtered to breast cancer cell lines via `primary_disease == "Breast Cancer"`
- Expression column: `rna_expression` (TPM)

**Results:**
- Breast cancer cell lines identified: confirmed panel
- Validated genes found in depmap: confirmed from 502

**Candidate Gene Expression Summary:**
| Gene | Mean TPM | Cell Line Tier |
|------|----------|----------------|
| STMN1 | ~7.5 | High |
| KPNA2 | ~7.0 | High |
| MAL2 | ~7.0 | High |
| TOP2A | ~6.5 | High ✅ Selected |
| CDK1 | ~6.2 | High ✅ Selected |
| PARP1 | ~6.5 | High |
| PLK1 | ~5.6 | Moderate ✅ Selected |
| SHCBP1 | ~3.8 | Moderate |
| ULBP2 | ~2.0 | Low |
| SPP1 | ~0.8 | Very Low |
| MMP13 | ~0.5 | Very Low |

> SPP1 and MMP13 showed near-zero cell line expression despite high patient tumor expression — indicating stromal/microenvironment origin rather than intrinsic tumor cell expression. These were deprioritized.

---

### `Phase7_Top30_CCLE_Expression.png`
> *Horizontal bar chart of top 30 validated genes by mean TPM in breast cancer cell lines. Color gradient from light orange (low) to dark red (high). Shows which of the 502 validated genes are most robustly expressed in experimental cell line models.*

###  `Phase7_Candidate_Genes_Boxplot.png`
> *Horizontal boxplot showing TPM distribution across breast cancer cell lines for 10 candidate genes. Each box represents variation across cell lines — wide boxes (ERBB2) indicate heterogeneous expression reflecting HER2-amplified vs non-amplified lines. Genes ordered by median expression. MMP13 and SPP1 visibly lowest — justifying their deprioritization.*

---

**Output Files:**
- `Phase7_CCLE_Mean_TPM_All_Validated.csv`
- `Phase7_CCLE_Candidate_Genes_TPM.csv`

---

## Phase 8 — Target Prioritization Scoring

**Objective:** Objectively rank all 502 validated genes using a multi-criteria evidence scoring system to select the top 3 drug targets.

### Scoring Criteria (Maximum 16 points)

| Criterion | Max Points | Rationale |
|-----------|-----------|-----------|
| Mean log2FC (P1+P2) | 3 | Expression strength across discovery datasets |
| TCGA log2FC | 3 | Validation strength in large cohort |
| Subtype specificity | 2 | Aggressive subtype (Basal/HER2) enrichment |
| Survival significance | 3 | Cox regression padj and HR direction |
| CCLE expression | 3 | Cell line tractability |
| Druggability bonus | +2 | Known oncogenic/druggable gene families |
| Stromal penalty | -3 | ECM/stromal genes deprioritized |

### Druggability Bonus Genes (selected list)
PARP1, ERBB2, KPNA2, STMN1, AURKB, PLK1, CDK1, BUB1, CCNB1, TOP2A, TYMS, RRM2, STAT1, FOXM1, KIF11, BIRC5, and others.

### Stromal Penalty Genes (selected list)
MMP1, MMP13, MMP14, COL1A1, COL1A2, COL3A1, COL5A1, COL5A2, FN1, POSTN, VCAN, BGN, SPP1, SULF1, SULF2, and others.

### Top 10 After Druggability-Adjusted Scoring
| Rank | Gene | Subtype | Mean log2FC | TCGA log2FC | CCLE TPM | Total Score |
|------|------|---------|-------------|-------------|----------|-------------|
| 1 | SHCBP1 | HER2-enriched | 3.09 | 2.64 | 3.81 | — |
| 2 | S100P | HER2-enriched | 4.66 | 4.82 | 6.28 | — |
| 3 | MYBL2 | Basal-like | 3.75 | 3.85 | 6.04 | — |
| 4 | FOXM1 | Basal-like | 3.10 | 3.58 | 5.65 | — |
| 5 | **PLK1** | Basal-like | 3.29 | 3.56 | 5.63 | — |
| 6 | KIF2C | Basal-like | 3.04 | 3.48 | 5.16 | — |
| 7 | KIF20A | Basal-like | 3.38 | 3.47 | 5.18 | — |
| 8 | **TOP2A** | HER2-enriched | 3.67 | 3.47 | 6.56 | — |
| 9 | RRM2 | HER2-enriched | 3.27 | 3.35 | 6.41 | — |
| 10 | **CDK1** | Basal-like | 3.02 | 3.16 | 6.24 | — |

> FOXM1 and MYBL2 excluded from final selection despite high scores — both are transcription factors with poorly defined binding pockets, making them unsuitable for structure-based molecular docking.

---

### `Phase8_Target_Prioritization_Scores.png`
> *Horizontal bar chart of top 20 scored genes. Bars color-coded by top PAM50 subtype (red=Basal-like, orange=HER2-enriched, blue=Luminal B). Score labels shown at end of each bar. PLK1, TOP2A, CDK1 visible in upper portion of chart.*

### `Phase8_Score_Breakdown_Heatmap.png`
> *Heatmap showing score breakdown for top 20 genes across 5 criteria columns (Fold Change, TCGA Validation, Subtype Specificity, Survival, Cell Line). Numbers displayed in each cell. Color from white (0) to dark red (max). Allows visual comparison of which evidence dimension drives each gene's total score.*

---

**Output Files:**
- `Phase8_Target_Prioritization_Scores.csv`
- `Phase8_Target_Prioritization_Scores_v2.csv` — druggability-adjusted final version

---

## Phase 9 — Druggability & Protein Structure Assessment

**Objective:** Compile complete structural biology profiles for the three selected targets to enable molecular docking setup.

---

### Target 1 — PLK1 (Polo-like Kinase 1)

| Attribute | Detail |
|-----------|--------|
| UniProt ID | P53350 |
| Gene | PLK1 |
| Function | Serine/threonine kinase; regulates centrosome maturation, spindle assembly, cytokinesis |
| Cancer Relevance | Overexpressed in TNBC/Basal-like; drives uncontrolled mitotic entry |
| AlphaFold ID | AF-P53350-F1 |

**PDB Structures:**
| PDB ID | Resolution | Description | Use |
|--------|-----------|-------------|-----|
| 2RKU | 2.1 Å | Kinase domain + ATP analog | ATP site docking |
| 3HIH | 1.9 Å | PLK1 + volasertib | Reference compound validation |
| 4O9W | 1.8 Å | Polo-box domain | PBD site docking |

**Active Site Residues:** Lys82, Asp176, Phe183, Cys67, Thr210

**Docking Box:** Center on Cys67–Lys82 region, size 20×20×20 Å

**Known Inhibitors:**
| Inhibitor | ChEMBL ID | Clinical Status |
|-----------|-----------|-----------------|
| Volasertib | CHEMBL1946170 | Phase III |
| BI-2536 | CHEMBL364881 | Phase II |
| Rigosertib | CHEMBL1614701 | Clinical |

---

### Target 2 — TOP2A (DNA Topoisomerase II Alpha)

| Attribute | Detail |
|-----------|--------|
| UniProt ID | P11388 |
| Gene | TOP2A |
| Function | Controls DNA topology; creates transient double-strand breaks during replication and transcription |
| Cancer Relevance | Overexpressed in HER2-enriched tumors; targeted by anthracyclines in clinical practice |
| AlphaFold ID | AF-P11388-F1 |

**PDB Structures:**
| PDB ID | Resolution | Description | Use |
|--------|-----------|-------------|-----|
| 4FM9 | 2.2 Å | Human TOP2A + etoposide | ✅ Primary docking structure |
| 1ZXM | 2.5 Å | ATPase domain | ATP site docking |
| 5GWK | 2.8 Å | TOP2A + doxorubicin | Intercalation site reference |

**Active Site Residues:** Tyr805 (cleavage site), Arg487, Lys489, Met762, Gln778

**Docking Box:** DNA cleavage gate around Tyr805, size 25×25×25 Å

**Known Inhibitors:**
| Inhibitor | ChEMBL ID | Clinical Status |
|-----------|-----------|-----------------|
| Etoposide | CHEMBL44657 | FDA Approved |
| Doxorubicin | CHEMBL53463 | FDA Approved |
| Dexrazoxane | CHEMBL1536 | FDA Approved |

---

### Target 3 — CDK1 (Cyclin-Dependent Kinase 1)

| Attribute | Detail |
|-----------|--------|
| UniProt ID | P06493 |
| Gene | CDK1 |
| Function | Master G2/M transition kinase; forms complex with Cyclin B1; essential for mitotic entry |
| Cancer Relevance | Cannot be bypassed — Basal-like/TNBC cells are particularly dependent due to RB/p53 loss |
| AlphaFold ID | AF-P06493-F1 |

**PDB Structures:**
| PDB ID | Resolution | Description | Use |
|--------|-----------|-------------|-----|
| 4YC6 | 2.0 Å | CDK1/Cyclin B complex | ✅ Physiological complex — primary structure |
| 6GU6 | 1.9 Å | CDK1 + RO-3306 | ✅ Selective inhibitor reference |
| 4Y72 | 2.1 Å | CDK1 + purvalanol | Alternative validation |

**Active Site Residues:** Lys33, Asp127, Phe80, Glu51, Leu83, Ala31

**Docking Box:** ATP binding cleft between N and C lobes around Lys33–Asp127 hinge, size 20×20×20 Å

**Known Inhibitors:**
| Inhibitor | ChEMBL ID | Clinical Status |
|-----------|-----------|-----------------|
| RO-3306 | CHEMBL1629769 | Research tool (selective CDK1) |
| Dinaciclib | CHEMBL1789844 | Phase III |
| Flavopiridol | CHEMBL13 | Clinical |

---

### Recommended Docking Workflow

**Software:** AutoDock Vina (open source, widely published)

1. Download PDB structures from [rcsb.org](https://www.rcsb.org) using IDs above
2. Prepare protein in AutoDockTools: remove waters, add polar hydrogens, assign Gasteiger charges, save as PDBQT
3. Download reference inhibitor SDF from [ChEMBL](https://www.ebi.ac.uk/chembl/), convert to PDBQT via Open Babel
4. Define grid box centered on active site residues listed above
5. Run docking: `exhaustiveness = 16` for publication quality
6. **Validate setup:** Redock co-crystallized ligand — RMSD must be < 2.0 Å before screening
7. Screen novel compounds against all three targets using validated setup
8. Analyze: binding affinity (kcal/mol), key interactions (H-bonds, hydrophobic contacts), RMSD

---

## Final Selected Targets

| Target | UniProt | PDB (Primary) | Subtype | log2FC | CCLE TPM | Key Advantage |
|--------|---------|---------------|---------|--------|----------|---------------|
| **PLK1** | P53350 | 3HIH | Basal-like | 3.29 | 5.63 | Deep ATP pocket, Phase III inhibitors as reference |
| **TOP2A** | P11388 | 4FM9 | HER2-enriched | 3.67 | 6.56 | FDA-approved inhibitors, highest CCLE expression |
| **CDK1** | P06493 | 6GU6 | Basal-like | 3.02 | 6.24 | Essential mitotic gate, cannot be bypassed by tumor cells |

---

## Output Files

### Complete File List by Phase

```
Project 3/
├── Phase1_Project1_Upregulated_Clean.csv
├── Phase1_Project2_Upregulated_Clean.csv
├── Phase2_Venn_Diagram_Overlap.png
├── Overlap_Genes.csv
├── Phase3_GO_BP_Dotplot.png
├── Phase3_GO_MF_Dotplot.png
├── Phase3_KEGG_Barplot.png
├── Phase3_Reactome_Dotplot.png
├── Phase3_GO_BP_Emapplot.png
├── Phase3_GO_BP_Results.csv
├── Phase3_GO_MF_Results.csv
├── Phase3_GO_CC_Results.csv
├── Phase3_KEGG_Results.csv
├── Phase3_Reactome_Results.csv
├── Phase4_TCGA_log2FC_Distribution.png
├── Phase4_Top30_Validated_Genes.png
├── Phase4_Tumor_vs_Normal_Scatter.png
├── Phase4_All_Genes_TCGA_Validation.csv
├── Phase4_Overlap_TCGA_Validated_Full.csv
├── Phase4_Confirmed_Upregulated_Genes.csv
├── Phase5_Subtype_Distribution.png
├── Phase5_Top20_Subtype_Boxplots.png
├── Phase5_Subtype_Heatmap.png
├── Phase5_Subtype_Mean_Expression.csv
├── Phase5_Gene_Top_Subtype.csv
├── Phase5_Subtype_Enrichment_Count.csv
├── Phase6_KM_PLK1.png
├── Phase6_KM_TOP2A.png
├── Phase6_KM_CDK1.png
├── Phase6_KM_PARP1.png
├── Phase6_KM_ERBB2.png
├── Phase6_KM_KPNA2.png
├── Phase6_KM_STAT1.png
├── Phase6_KM_SPP1.png
├── Phase6_Forest_Plot.png
├── Phase6_Cox_Survival_Results.csv
├── Phase7_Top30_CCLE_Expression.png
├── Phase7_Candidate_Genes_Boxplot.png
├── Phase7_CCLE_Mean_TPM_All_Validated.csv
├── Phase7_CCLE_Candidate_Genes_TPM.csv
├── Phase8_Target_Prioritization_Scores.csv
├── Phase8_Target_Prioritization_Scores_v2.csv
└── Phase8_Score_Breakdown_Heatmap.png
```

---

## Key Results Summary

| Phase | Key Finding |
|-------|-------------|
| Phase 1–2 | 666 genes consistently upregulated across two independent breast cancer datasets |
| Phase 3 | Cell cycle dysregulation and chromosomal instability are dominant biological themes |
| Phase 4 | 93% validation rate in full TCGA-BRCA cohort (502/540 genes confirmed) |
| Phase 5 | 502-gene signature enriched in Basal-like (44%) and HER2-enriched (22%) subtypes |
| Phase 6 | 4 genes survive Cox multiple testing correction; low event rate limits power |
| Phase 7 | MMP13 and SPP1 deprioritized — near-zero cell line expression indicates stromal origin |
| Phase 8 | Druggability-adjusted scoring surfaces enzymatic tumor-intrinsic targets |
| Phase 9 | PLK1, TOP2A, CDK1 profiled with PDB structures, active sites, reference inhibitors |

---

## Citation & Acknowledgements

Data sources:
- TCGA Research Network: https://www.cancer.gov/tcga
- curatedTCGAData: Ramos et al., Cancer Research (2020)
- DepMap: Broad Institute Cancer Dependency Map
- Reactome: Fabregat et al., Nucleic Acids Research (2018)
- KEGG: Kanehisa & Goto, Nucleic Acids Research (2000)

---

*Pipeline developed in R 4.5.2 | Bioconductor 3.22 | Generated as part of PhD research in computational oncology*
