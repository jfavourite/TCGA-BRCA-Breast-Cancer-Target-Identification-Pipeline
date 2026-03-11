# ============================================================
# Project 3: Intersection, Validation & Target Prioritization
# TCGA-BRCA Breast Cancer
# ============================================================
install.packages("ggVennDiagram")
library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(ReactomePA)
library(ggplot2)
library(forcats)
library(dplyr)
library(readr)
library(ggplot2)
library(ggVennDiagram)
library(curatedTCGAData)
library(MultiAssayExperiment)
library(TCGAutils)
library(ggplot2)
library(ggpubr)
library(tidyr)
library(ggpubr)
library(pheatmap)
library(survival)
library(survminer)
library(depmap)
library(ExperimentHub)

#  Load both CSV files from previous projects

p1_raw <- read_csv("TCGA_BRCA_HighConfidence_DEGs1_Project_1.csv")
p2_raw <- read_csv("WGS_Breast_cancer_Full_Differential_Expression_project_2.csv")

p1_raw <- p1_raw %>% rename(symbol = Symbol)


# Removes the .version number e.g. ENSG00000163092.21 → ENSG00000163092
p1_raw <- p1_raw %>%
  mutate(ENSG_base = sub("\\..*", "", ENSG))

p2_raw <- p2_raw %>%
  mutate(ENSG_base = sub("\\..*", "", ENSG))

p1_na <- p1_raw %>%
  filter(is.na(log2FoldChange) | is.na(padj) | is.na(symbol))

p2_na <- p2_raw %>%
  filter(is.na(log2FoldChange) | is.na(padj) | is.na(symbol))

p1_clean <- p1_raw %>%
  filter(!is.na(log2FoldChange), !is.na(padj), !is.na(symbol))

p2_clean <- p2_raw %>%
  filter(!is.na(log2FoldChange), !is.na(padj), !is.na(symbol))

# Apply DEG thresholds


p1_up <- p1_clean %>%
  filter(log2FoldChange > 1, padj < 0.05)

p2_up <- p2_clean %>%
  filter(log2FoldChange > 1, padj < 0.05)

#  Remove duplicates

cat("Project 1 — duplicate symbols before removal:",
    sum(duplicated(p1_up$symbol)), "\n")
cat("Project 2 — duplicate symbols before removal:",
    sum(duplicated(p2_up$symbol)), "\n\n")

# Keep the row with lowest padj per symbol (most significant)
p1_up <- p1_up %>%
  arrange(padj) %>%
  distinct(symbol, .keep_all = TRUE)

p2_up <- p2_up %>%
  arrange(padj) %>%
  distinct(symbol, .keep_all = TRUE)

#  Extract final clean gene vectors
project1_up <- p1_up$symbol
project2_up <- p2_up$symbol


#  Save clean files
write_csv(p1_up, "Phase1_Project1_Upregulated_Clean.csv")
write_csv(p2_up, "Phase1_Project2_Upregulated_Clean.csv")

overlap_genes <- intersect(project1_up, project2_up)

print(overlap_genes)

 # Venn Diagram
venn_list <- list(
  "Project 1\n(TCGA-BRCA\nTumor vs Normal)" = project1_up,
  "Project 2\n(Salmon\nQuantification)"      = project2_up
)

venn_plot <- ggVennDiagram(
  venn_list,
  label_alpha = 0,           # transparent label background
  label = "count"            # show counts in each region
) +
  scale_fill_gradient(
    low  = "#d4eaf7",
    high = "#1a6fa8"
  ) +
  scale_color_manual(values = c("#1a6fa8", "#e07b39")) +
  labs(
    title    = "Upregulated Gene Overlap — Project 1 vs Project 2",
    subtitle = "TCGA-BRCA Breast Cancer | log2FC > 1, padj < 0.05",
    caption  = paste0("Intersection: ", length(overlap_genes), " high-confidence tumor-upregulated genes")
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 14, hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40"),
    plot.caption  = element_text(size = 10, hjust = 0.5, color = "grey30"),
    legend.title  = element_text(size = 9),
    legend.position = "right"
  )

print(venn_plot)

# Save Venn diagram
ggsave(
  "Phase2_Venn_Diagram_Overlap.png",
  plot   = venn_plot,
  width  = 8,
  height = 6,
  dpi    = 300,
  bg     = "white"
)


# Build Overlap Gene Table

p1_overlap <- p1_up %>%
  filter(symbol %in% overlap_genes) %>%
  select(symbol, ENSG_base,
         log2FC_P1  = log2FoldChange,
         padj_P1    = padj)

p2_overlap <- p2_up %>%
  filter(symbol %in% overlap_genes) %>%
  select(symbol,
         log2FC_P2  = log2FoldChange,
         padj_P2    = padj)

overlap_table <- p1_overlap %>%
  left_join(p2_overlap, by = "symbol") %>%
  mutate(
    mean_log2FC = round((log2FC_P1 + log2FC_P2) / 2, 3),
    min_padj    = pmin(padj_P1, padj_P2)
  ) %>%
  arrange(desc(mean_log2FC))

print(head(overlap_table, 20))

#  Save Overlap CSV
write_csv(overlap_table, "Overlap_Genes.csv")

overlap_symbols <- overlap_table$symbol


#  — Functional Characterization of Overlap Genes
# GO, KEGG & Reactome Enrichment Analysis
#  Convert Gene Symbols → Entrez IDs 

symbol_to_entrez <- bitr(
  overlap_symbols,
  fromType = "SYMBOL",
  toType   = "ENTREZID",
  OrgDb    = org.Hs.eg.db
)

entrez_ids <- symbol_to_entrez$ENTREZID

#. GO Enrichment Analysis 

# Biological Process
go_bp <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE   
)

# Molecular Function
go_mf <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# Cellular Component
go_cc <- enrichGO(
  gene          = entrez_ids,
  OrgDb         = org.Hs.eg.db,
  ont           = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

#  KEGG Enrichment Analysis

kegg_res <- enrichKEGG(
  gene          = entrez_ids,
  organism      = "hsa",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05
)

# Convert Entrez to symbols for readability
kegg_res <- setReadable(kegg_res, OrgDb = org.Hs.eg.db,
                        keyType = "ENTREZID")

# Reactome Enrichment Analysis

reactome_res <- enrichPathway(
  gene          = entrez_ids,
  organism      = "human",
  pvalueCutoff  = 0.05,
  qvalueCutoff  = 0.05,
  readable      = TRUE
)

# VISUALIZATIONS

# a. GO Biological Process — Dotplot
bp_dot <- dotplot(go_bp, showCategory = 20) +
  labs(
    title    = "GO Biological Process Enrichment",
    subtitle = "Top 20 terms | 666 Overlap Genes | TCGA-BRCA",
    x        = "Gene Ratio"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y   = element_text(size = 8)
  )

print(bp_dot)
ggsave("Phase3_GO_BP_Dotplot.png", bp_dot,
       width = 10, height = 8, dpi = 300, bg = "white")

# b. GO Molecular Function — Dotplot
mf_dot <- dotplot(go_mf, showCategory = 15) +
  labs(
    title    = "GO Molecular Function Enrichment",
    subtitle = "Top 15 terms | 666 Overlap Genes | TCGA-BRCA",
    x        = "Gene Ratio"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y   = element_text(size = 8)
  )

print(mf_dot)
ggsave("Phase3_GO_MF_Dotplot.png", mf_dot,
       width = 10, height = 7, dpi = 300, bg = "white")

# c. KEGG — Barplot
kegg_bar <- barplot(kegg_res, showCategory = 20) +
  labs(
    title    = "KEGG Pathway Enrichment",
    subtitle = "Top 20 pathways | 666 Overlap Genes | TCGA-BRCA",
    x        = "Gene Count"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y   = element_text(size = 8)
  )

print(kegg_bar)
ggsave("Phase3_KEGG_Barplot.png", kegg_bar,
       width = 10, height = 8, dpi = 300, bg = "white")

# d. Reactome — Dotplot
reactome_dot <- dotplot(reactome_res, showCategory = 20) +
  labs(
    title    = "Reactome Pathway Enrichment",
    subtitle = "Top 20 pathways | 666 Overlap Genes | TCGA-BRCA",
    x        = "Gene Ratio"
  ) +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y   = element_text(size = 8)
  )

print(reactome_dot)
ggsave("Phase3_Reactome_Dotplot.png", reactome_dot,
       width = 10, height = 8, dpi = 300, bg = "white")

# e. GO-BP Enrichment Map (term relationships)
go_bp2 <- pairwise_termsim(go_bp)
emap <- emapplot(go_bp2, showCategory = 30) +
  labs(title = "GO-BP Enrichment Map — Term Relationships") +
  theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5))

print(emap)
ggsave("Phase3_GO_BP_Emapplot.png", emap,
       width = 12, height = 10, dpi = 300, bg = "white")

#  Save Result Tables
write_csv(as.data.frame(go_bp),      "Phase3_GO_BP_Results.csv")
write_csv(as.data.frame(go_mf),      "Phase3_GO_MF_Results.csv")
write_csv(as.data.frame(go_cc),      "Phase3_GO_CC_Results.csv")
write_csv(as.data.frame(kegg_res),   "Phase3_KEGG_Results.csv")
write_csv(as.data.frame(reactome_res),"Phase3_Reactome_Results.csv")

print(go_bp@result %>%
        filter(p.adjust < 0.05) %>%
        dplyr::select(Description, GeneRatio, p.adjust, Count) %>%
        head(10))

print(kegg_res@result %>%
        filter(p.adjust < 0.05) %>%
        dplyr::select(Description, GeneRatio, p.adjust, Count) %>%
        head(10))
print(reactome_res@result %>%
        filter(p.adjust < 0.05) %>%
        dplyr::select(Description, GeneRatio, p.adjust, Count) %>%
        head(10))


# — TCGA Validation via curatedTCGAData
# External validation of 666 overlap genes in TCGA-BRCA

#  Stream TCGA-BRCA RNASeq2 Data

brca_mae <- curatedTCGAData(
  diseaseCode = "BRCA",
  assays      = "RNASeq2GeneNorm",
  version     = "2.0.1",
  dry.run     = FALSE
)

print(brca_mae)

# Extract Expression Matrix
# Extract the RNASeq2GeneNorm assay
rna_assay <- assay(brca_mae[[1]])

# Separate Tumor vs Normal Samples
# TCGA barcode positions 14-15 indicate sample type:
# 01 = Primary Solid Tumor
# 11 = Solid Tissue Normal
sample_ids   <- colnames(rna_assay)
sample_types <- substr(sample_ids, 14, 15)

tumor_ids  <- sample_ids[sample_types == "01"]
normal_ids <- sample_ids[sample_types == "11"]

# Subset matrix
tumor_mat  <- rna_assay[, tumor_ids]
normal_mat <- rna_assay[, normal_ids]

# Filter to Overlap Genes Only
genes_in_tcga <- intersect(overlap_symbols, rownames(rna_assay))

# Subset to overlap genes
tumor_sub  <- tumor_mat[genes_in_tcga, ]
normal_sub <- normal_mat[genes_in_tcga, ]

# Compute Mean Expression Per Gene and comarison table

mean_tumor  <- rowMeans(tumor_sub,  na.rm = TRUE)
mean_normal <- rowMeans(normal_sub, na.rm = TRUE)

validation_df <- data.frame(
  symbol      = genes_in_tcga,
  mean_tumor  = round(mean_tumor,  3),
  mean_normal = round(mean_normal, 3),
  log2FC_tcga = round(log2((mean_tumor + 1) / (mean_normal + 1)), 3)
) %>%
  arrange(desc(log2FC_tcga))

# Confirm Overexpression in Tumor

validated_up <- validation_df %>%
  filter(log2FC_tcga > 1)

validated_not <- validation_df %>%
  filter(log2FC_tcga <= 1)

val_rate <- round(nrow(validated_up) / length(genes_in_tcga) * 100, 1)

#  Merge with Original Overlap Table
overlap_validated <- overlap_table %>%
  filter(symbol %in% genes_in_tcga) %>%
  left_join(validation_df, by = "symbol") %>%
  arrange(desc(log2FC_tcga))

# Flag validation status
overlap_validated <- overlap_validated %>%
  mutate(tcga_validated = ifelse(log2FC_tcga > 1, "YES", "NO"))

print(overlap_validated %>%
        dplyr::select(symbol, log2FC_P1, log2FC_P2,
                      log2FC_tcga, tcga_validated) %>%
        head(20))

#  Visualizations

#  log2FC Distribution — Tumor vs Normal
fc_plot <- ggplot(validation_df, aes(x = log2FC_tcga)) +
  geom_histogram(binwidth = 0.2, fill = "#1a6fa8",
                 color = "white", alpha = 0.85) +
  geom_vline(xintercept = 1, linetype = "dashed",
             color = "red", linewidth = 0.8) +
  annotate("text", x = 1.3, y = Inf, vjust = 2,
           label = "log2FC = 1 threshold",
           color = "red", size = 3.5) +
  labs(
    title    = "TCGA-BRCA: log2FC Distribution of Overlap Genes",
    subtitle = "Tumor vs Normal | 666 Overlap Genes",
    x        = "log2 Fold Change (Tumor / Normal)",
    y        = "Number of Genes"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40")
  )

print(fc_plot)
ggsave("Phase4_TCGA_log2FC_Distribution.png", fc_plot,
       width = 9, height = 6, dpi = 300, bg = "white")

#  Top 30 Validated Genes — Barplot
top30 <- overlap_validated %>%
  filter(tcga_validated == "YES") %>%
  head(30)

bar_plot <- ggplot(top30,
                   aes(x = reorder(symbol, log2FC_tcga),
                       y = log2FC_tcga,
                       fill = log2FC_tcga)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = "#74b9e8", high = "#1a3a6b") +
  coord_flip() +
  labs(
    title    = "Top 30 TCGA-Validated Upregulated Genes",
    subtitle = "TCGA-BRCA Tumor vs Normal | log2FC > 1",
    x        = NULL,
    y        = "log2 Fold Change (TCGA)",
    fill     = "log2FC"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.y   = element_text(size = 8)
  )

print(bar_plot)
ggsave("Phase4_Top30_Validated_Genes.png", bar_plot,
       width = 9, height = 8, dpi = 300, bg = "white")

# Scatterplot — mean tumor vs mean normal expression
scatter_plot <- ggplot(validation_df,
                       aes(x = log2(mean_normal + 1),
                           y = log2(mean_tumor  + 1),
                           color = log2FC_tcga > 1)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  scale_color_manual(
    values = c("FALSE" = "grey70", "TRUE" = "#1a6fa8"),
    labels = c("FALSE" = "log2FC ≤ 1", "TRUE" = "log2FC > 1 (validated)")
  ) +
  labs(
    title    = "Tumor vs Normal Expression — TCGA-BRCA",
    subtitle = "Each point = one overlap gene",
    x        = "log2 Mean Normal Expression",
    y        = "log2 Mean Tumor Expression",
    color    = "Validation Status"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "bottom"
  )

print(scatter_plot)
ggsave("Phase4_Tumor_vs_Normal_Scatter.png", scatter_plot,
       width = 8, height = 7, dpi = 300, bg = "white")

# Save Outputs
write_csv(validation_df,       "Phase4_All_Genes_TCGA_Validation.csv")
write_csv(overlap_validated,   "Phase4_Overlap_TCGA_Validated_Full.csv")
write_csv(validated_up,        "Phase4_Confirmed_Upregulated_Genes.csv")


validated_symbols <- overlap_validated %>%
  filter(tcga_validated == "YES") %>%
  pull(symbol)

# — Subtype Analysis
# Luminal A, Luminal B, HER2-enriched, Basal-like/TNBC

# Extract Clinical & Subtype Data

clinical <- colData(brca_mae) %>% as.data.frame()
print(names(clinical))

# ── 2. Find PAM50 Subtype Column
pam50_col <- names(clinical)[grep("PAM50|subtype|Subtype|intrinsic",
                                  names(clinical),
                                  ignore.case = TRUE)]
print(pam50_col)

# ── 3. Extract & Clean Subtype Labels ────────────────────────
# Use the first PAM50 column found
subtype_col <- pam50_col[1]

subtype_df <- clinical %>%
  dplyr::select(patientID, subtype = all_of(subtype_col)) %>%
  filter(!is.na(subtype)) %>%
  mutate(subtype = gsub("BRCA\\.|Subtype\\.", "", subtype)) 

print(table(subtype_df$subtype))

# Match Samples Between Expression & Subtype

# Get tumor expression for validated genes only
tumor_validated <- tumor_mat[validated_symbols, ]

# TCGA barcodes in expression matrix need to be trimmed to
# match patient IDs in clinical data (first 12 characters)
colnames_trimmed <- substr(colnames(tumor_validated), 1, 12)
colnames(tumor_validated) <- colnames_trimmed

# Keep only samples that have subtype annotation
common_samples <- intersect(colnames(tumor_validated),
                            subtype_df$patientID)

# Subset expression and subtype to common samples
expr_sub  <- tumor_validated[, common_samples]
subtype_matched <- subtype_df %>%
  filter(patientID %in% common_samples) %>%
  arrange(match(patientID, common_samples))

# Compute Mean Expression Per Gene Per Subtype
# Build a long-format table: gene × sample → add subtype
expr_long <- as.data.frame(t(expr_sub)) %>%
  mutate(patientID = rownames(.)) %>%
  left_join(subtype_df, by = "patientID") %>%
  filter(!is.na(subtype)) %>%
  pivot_longer(cols      = -c(patientID, subtype),
               names_to  = "symbol",
               values_to = "expression")

# Mean expression per gene per subtype
subtype_means <- expr_long %>%
  group_by(symbol, subtype) %>%
  summarise(mean_expr = mean(expression, na.rm = TRUE),
            .groups = "drop")

print(head(top_subtype_per_gene))
# For each gene, find which subtype has highest mean expression
top_subtype_per_gene <- subtype_means %>%
  group_by(symbol) %>%
  dplyr::slice_max(mean_expr, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::rename(top_subtype = subtype, top_mean_expr = mean_expr) %>%
  as.data.frame()
print(head(top_subtype_per_gene))

# Count how many genes are enriched per subtype
subtype_enrichment_count <- top_subtype_per_gene %>%
  dplyr::count(top_subtype, name = "n_genes") %>%
  arrange(desc(n_genes))

print(subtype_enrichment_count)

# Visualizations

# a. Subtype sample distribution pie/bar
subtype_bar <- ggplot(subtype_df %>%
                        filter(patientID %in% common_samples),
                      aes(x = subtype, fill = subtype)) +
  geom_bar(color = "white") +
  scale_fill_manual(values = c(
    "LumA"    = "#2196F3",
    "LumB"    = "#4CAF50",
    "Her2"    = "#FF9800",
    "Basal"   = "#F44336",
    "Normal"  = "#9C27B0"
  )) +
  labs(
    title    = "TCGA-BRCA Sample Distribution by PAM50 Subtype",
    subtitle = "Samples with matched subtype annotation",
    x        = "PAM50 Subtype",
    y        = "Number of Samples",
    fill     = "Subtype"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    legend.position = "none"
  )

print(subtype_bar)
ggsave("Phase5_Subtype_Distribution.png", subtype_bar,
       width = 8, height = 6, dpi = 300, bg = "white")

# b. Top 20 genes — expression across subtypes (boxplot)
# Pick top 20 genes by overall mean expression
top20_genes <- subtype_means %>%
  group_by(symbol) %>%
  summarise(overall_mean = mean(mean_expr)) %>%
  slice_max(overall_mean, n = 20) %>%
  pull(symbol)

box_data <- expr_long %>%
  filter(symbol %in% top20_genes)

box_plot <- ggplot(box_data,
                   aes(x = subtype, y = log2(expression + 1),
                       fill = subtype)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.85) +
  facet_wrap(~symbol, scales = "free_y", ncol = 4) +
  scale_fill_manual(values = c(
    "LumA"   = "#2196F3",
    "LumB"   = "#4CAF50",
    "Her2"   = "#FF9800",
    "Basal"  = "#F44336",
    "Normal" = "#9C27B0"
  )) +
  labs(
    title    = "Top 20 Validated Genes — Expression Across PAM50 Subtypes",
    subtitle = "TCGA-BRCA | log2(expression + 1)",
    x        = NULL,
    y        = "log2(Expression + 1)",
    fill     = "Subtype"
  ) +
  theme_classic() +
  theme(
    plot.title    = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x   = element_text(angle = 45, hjust = 1, size = 7),
    strip.text    = element_text(face = "bold", size = 8),
    legend.position = "bottom"
  )

print(box_plot)
ggsave("Phase5_Top20_Subtype_Boxplots.png", box_plot,
       width = 14, height = 12, dpi = 300, bg = "white")

# c. Heatmap — top 40 genes × subtypes
top40_genes <- subtype_means %>%
  group_by(symbol) %>%
  summarise(overall_mean = mean(mean_expr)) %>%
  slice_max(overall_mean, n = 40) %>%
  pull(symbol)

heatmap_mat <- subtype_means %>%
  filter(symbol %in% top40_genes) %>%
  pivot_wider(names_from  = subtype,
              values_from = mean_expr) %>%
  tibble::column_to_rownames("symbol") %>%
  as.matrix()

# Scale rows (z-score per gene)
heatmap_scaled <- t(scale(t(heatmap_mat)))

png("Phase5_Subtype_Heatmap.png",
    width = 10, height = 12,
    units = "in", res = 300)

pheatmap(
  heatmap_scaled,
  color            = colorRampPalette(c("#2166ac","white","#d6604d"))(100),
  cluster_cols     = TRUE,
  cluster_rows     = TRUE,
  show_rownames    = TRUE,
  fontsize_row     = 7,
  fontsize_col     = 10,
  main             = "Top 40 Validated Genes — Mean Expression by PAM50 Subtype\n(z-score scaled)",
  border_color     = NA,
  annotation_legend = TRUE
)

# Save Outputs
write_csv(subtype_means,          "Phase5_Subtype_Mean_Expression.csv")
write_csv(top_subtype_per_gene,   "Phase5_Gene_Top_Subtype.csv")
write_csv(subtype_enrichment_count, "Phase5_Subtype_Enrichment_Count.csv")

# — Survival Analysis
subtype_for_survival <- subtype_matched

# Extract survival data from TCGA
clinical_survival <- as.data.frame(colData(brca_mae))

print(names(clinical_survival)[grepl("vital|death|follow|days|status|surv", 
                                     names(clinical_survival), ignore.case = TRUE)])

surv_cols <- names(clinical_survival)[grepl("vital|death|follow|days", 
                                            names(clinical_survival), ignore.case = TRUE)]
print(surv_cols)

# Build survival data frame
surv_df <- clinical_survival %>%
  dplyr::select(patientID, vital_status, days_to_death, days_to_last_followup) %>%
  mutate(
    os_time = ifelse(!is.na(days_to_death) & days_to_death > 0, 
                     days_to_death, days_to_last_followup),
    os_event = as.numeric(as.character(vital_status))  # converts "1" → 1, "0" → 0
  ) %>%
  filter(!is.na(os_time), os_time > 0, !is.na(os_event))

# Merge with expression data
# Use top candidates from heatmap
target_genes <- c("PARP1", "ERBB2", "KPNA2", "STAT1", "STMN1", "TUBA1C", "SPP1", "CRABP2")

expr_targets <- as.data.frame(t(tumor_validated[target_genes, ])) %>%
  mutate(patientID = substr(rownames(.), 1, 12))

surv_expr <- surv_df %>%
  inner_join(expr_targets, by = "patientID")

# KM curves per gene
km_plots <- list()

for (gene in target_genes) {
  
  df <- surv_expr %>%
    dplyr::select(os_time, os_event, expr = all_of(gene)) %>%
    filter(!is.na(expr)) %>%
    mutate(group = ifelse(expr >= median(expr), "High", "Low"))
  
  fit <- survfit(Surv(os_time, os_event) ~ group, data = df)
  
  p <- ggsurvplot(fit, data = df,
                  title = paste0(gene, " Expression — Overall Survival"),
                  legend.title = "Expression",
                  legend.labs = c("High", "Low"),
                  palette = c("#E41A1C", "#377EB8"),
                  pval = TRUE, pval.method = TRUE,
                  risk.table = TRUE,
                  xlab = "Days", ylab = "Survival Probability",
                  ggtheme = theme_bw(base_size = 13))
  
  km_plots[[gene]] <- p
}

# Save individual KM plots
for (gene in target_genes) {
  png(paste0("Phase6_KM_", gene, ".png"), width = 900, height = 700)
  print(km_plots[[gene]]$plot)
  dev.off()
}

# Cox regression for all 502 validated genes
cox_results <- list()

all_genes <- rownames(tumor_validated)
all_genes <- all_genes[all_genes %in% rownames(tumor_validated)]

expr_all <- as.data.frame(t(tumor_validated)) %>%
  mutate(patientID = substr(rownames(.), 1, 12))

surv_all <- surv_df %>%
  inner_join(expr_all, by = "patientID")

for (gene in validated_symbols) {
  if (!gene %in% names(surv_all)) next
  
  df <- surv_all %>%
    dplyr::select(os_time, os_event, expr = all_of(gene)) %>%
    filter(!is.na(expr))
  
  tryCatch({
    fit <- coxph(Surv(os_time, os_event) ~ expr, data = df)
    s <- summary(fit)
    cox_results[[gene]] <- data.frame(
      symbol    = gene,
      HR        = s$conf.int[1, 1],
      HR_lower  = s$conf.int[1, 3],
      HR_upper  = s$conf.int[1, 4],
      pvalue    = s$coefficients[1, 5],
      stringsAsFactors = FALSE
    )
  }, error = function(e) NULL)
}

cox_df <- bind_rows(cox_results) %>%
  mutate(padj = p.adjust(pvalue, method = "BH")) %>%
  arrange(pvalue)

print(head(cox_df, 10))

# Forest plot — top 15 significant genes
top_cox <- cox_df %>% filter(padj < 0.05) %>% head(15)

if (nrow(top_cox) > 0) {
  png("Phase6_Forest_Plot.png", width = 900, height = 700)
  print(
    ggplot(top_cox, aes(x = HR, y = reorder(symbol, -HR),
                        xmin = HR_lower, xmax = HR_upper)) +
      geom_pointrange(color = "#E41A1C", size = 0.7) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "gray40") +
      labs(title = "Cox Regression — Significant Survival Genes",
           x = "Hazard Ratio (HR)", y = NULL,
           caption = "HR > 1 = worse survival with high expression") +
      theme_bw(base_size = 13)
  )
  dev.off()
}
# Save results
write.csv(cox_df, "Phase6_Cox_Survival_Results.csv", row.names = FALSE)

# Cell Line Expression (depmap
# Load depmap data
eh <- ExperimentHub()
query(eh, "depmap")

# Load latest TPM expression + metadata
tpm_data <- eh[["EH7556"]]      # TPM_22Q2
meta_data <- eh[["EH7558"]]     # metadata_22Q2

print(names(meta_data))


# Filter to breast cancer cell lines
print(sort(unique(meta_data$primary_disease)))

breast_lines <- meta_data %>%
  filter(grepl("Breast", primary_disease, ignore.case = TRUE))

print(table(breast_lines$lineage_molecular_subtype))

# Filter TPM to breast lines + validated genes
breast_tpm <- tpm_data %>%
  filter(depmap_id %in% breast_lines$depmap_id) %>%
  filter(gene_name %in% validated_symbols)

#  Add metadata
breast_tpm <- breast_tpm %>%
  left_join(breast_lines %>% dplyr::select(depmap_id, cell_line_name, 
                                           lineage_molecular_subtype,
                                           primary_or_metastasis), 
            by = "depmap_id")

print(head(breast_tpm))
print(names(breast_tpm))

gene_mean_tpm <- breast_tpm %>%
  group_by(gene_name) %>%
  summarise(
    mean_tpm   = mean(rna_expression, na.rm = TRUE),
    median_tpm = stats::median(rna_expression, na.rm = TRUE),
    n_lines    = n_distinct(depmap_id),
    .groups    = "drop"
  ) %>%
  arrange(desc(mean_tpm))

print(head(gene_mean_tpm, 20))

#  Plots
top30_ccle <- head(gene_mean_tpm, 30)

png("Phase7_Top30_CCLE_Expression.png", width = 1000, height = 700)
ggplot(top30_ccle, aes(x = reorder(gene_name, mean_tpm), y = mean_tpm, fill = mean_tpm)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  scale_fill_gradient(low = "#fee8c8", high = "#b30000") +
  labs(title = "Top 30 Validated Genes — Mean TPM in Breast Cancer Cell Lines (CCLE)",
       x = NULL, y = "Mean TPM Expression", fill = "Mean TPM") +
  theme_bw(base_size = 13)

candidates <- c("PARP1", "ERBB2", "KPNA2", "STAT1", "STMN1", 
                "MMP13", "SHCBP1", "MAL2", "ULBP2", "SPP1")

candidate_tpm <- breast_tpm %>%
  filter(gene_name %in% candidates)

png("Phase7_Candidate_Genes_Boxplot.png", width = 1100, height = 700)
ggplot(candidate_tpm, aes(x = reorder(gene_name, rna_expression, stats::median), 
                          y = rna_expression, fill = gene_name)) +
  geom_boxplot(outlier.size = 0.8) +
  coord_flip() +
  labs(title = "Candidate Gene Expression Across Breast Cancer Cell Lines (CCLE)",
       x = NULL, y = "TPM Expression") +
  theme_bw(base_size = 13) +
  theme(legend.position = "none")

write.csv(gene_mean_tpm, "Phase7_CCLE_Mean_TPM_All_Validated.csv", row.names = FALSE)
write.csv(candidate_tpm, "Phase7_CCLE_Candidate_Genes_TPM.csv",    row.names = FALSE)



# Target Prioritization Scoring
# Build master evidence table

# Base: all 502 validated genes
scoring_base <- overlap_validated %>%
  filter(symbol %in% validated_symbols) %>%
  dplyr::select(symbol, mean_log2FC, min_padj)

# Add TCGA validation log2FC
tcga_fc <- read.csv("Phase4_Confirmed_Upregulated_Genes.csv") %>%
  dplyr::select(symbol, tcga_log2FC = log2FC_tcga)

# Add subtype info
subtype_info <- top_subtype_per_gene %>%
  dplyr::select(symbol, top_subtype, top_mean_expr)

# Add Cox survival
cox_info <- cox_df %>%
  dplyr::select(symbol, HR, cox_padj = padj)

# Add CCLE expression
ccle_info <- gene_mean_tpm %>%
  dplyr::select(symbol = gene_name, mean_tpm)

# ── Step 2: Merge all evidence ────────────────────────────────────────────────
master <- scoring_base %>%
  left_join(tcga_fc,     by = "symbol") %>%
  left_join(subtype_info, by = "symbol") %>%
  left_join(cox_info,    by = "symbol") %>%
  left_join(ccle_info,   by = "symbol")

cat("Master table dimensions:", dim(master), "\n")
cat("Columns:", names(master), "\n")



# ── Step 3: Build scoring system ──────────────────────────────────────────────

scored2 <- master %>%
  mutate(
    # 1. Expression strength (0-3 points)
    score_fc = case_when(
      mean_log2FC >= 4  ~ 3,
      mean_log2FC >= 3  ~ 2,
      mean_log2FC >= 2  ~ 1,
      TRUE              ~ 0
    ),
    
    # 2. TCGA validation strength (0-3 points)
    score_tcga = case_when(
      tcga_log2FC >= 4  ~ 3,
      tcga_log2FC >= 3  ~ 2,
      tcga_log2FC >= 2  ~ 1,
      TRUE              ~ 0
    ),
    
    # 3. Subtype specificity (0-2 points)
    score_subtype = case_when(
      top_subtype == "Basal-like"    ~ 2,
      top_subtype == "HER2-enriched" ~ 2,
      top_subtype == "Luminal B"     ~ 1,
      TRUE                           ~ 0
    ),
    
    # 4. Survival significance (0-3 points)
    score_survival = case_when(
      !is.na(cox_padj) & cox_padj < 0.05 & HR > 1 ~ 3,
      !is.na(cox_padj) & cox_padj < 0.05           ~ 2,
      !is.na(cox_padj) & cox_padj < 0.20           ~ 1,
      TRUE                                          ~ 0
    ),
    
    # 5. Cell line expression (0-3 points) — KEY FILTER
    score_ccle = case_when(
      mean_tpm >= 7  ~ 3,
      mean_tpm >= 5  ~ 2,
      mean_tpm >= 2  ~ 1,
      TRUE           ~ 0
    ),
    
    # 6. Druggability bonus (+2 points)
    score_druggable = ifelse(symbol %in% druggable_genes, 2, 0),
    
    # 7. Stromal penalty (-3 points)
    score_stromal = ifelse(symbol %in% stromal_genes, -3, 0),
    
    # ── Total score (max = 16) ────────────────────────────────────────────────
    total_score = score_fc + score_tcga + score_subtype + 
      score_survival + score_ccle + score_druggable + score_stromal
  ) %>%
  arrange(desc(total_score))

# ── Step 4: Visualization ─────────────────────────────────────────────────────
top20_scored <- head(scored2, 20)

png("Phase8_Target_Prioritization_Scores.png", width = 1000, height = 800)
ggplot(top20_scored, aes(x = reorder(symbol, total_score), y = total_score, fill = top_subtype)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = total_score), hjust = -0.3, size = 4) +
  coord_flip() +
  scale_fill_manual(values = c("Basal-like"    = "#E41A1C",
                               "HER2-enriched" = "#FF7F00",
                               "Luminal A"     = "#4DAF4A",
                               "Luminal B"     = "#377EB8",
                               "Normal-like"   = "#984EA3")) +
  labs(title = "Top 20 Prioritized Drug Targets — Evidence-Based Scoring",
       x = NULL, y = "Total Score (max = 14)", fill = "Top Subtype") +
  theme_bw(base_size = 13) +
  expand_limits(y = 16)
dev.off()

#  Score breakdown heatmap for top 20 ────────────────────────────────
score_mat <- top20_scored %>%
  dplyr::select(symbol, score_fc, score_tcga, score_subtype, 
                score_survival, score_ccle) %>%
  tibble::column_to_rownames("symbol") %>%
  as.matrix()

colnames(score_mat) <- c("Fold Change", "TCGA Validation", 
                         "Subtype Specificity", "Survival", "Cell Line")

pdf("Phase8_Score_Breakdown_Heatmap.pdf", width = 10, height = 12)
pheatmap(score_mat,
         cluster_rows = FALSE, cluster_cols = FALSE,
         display_numbers = TRUE, number_format = "%d",
         color = colorRampPalette(c("white", "#fdae61", "#d73027"))(10),
         main = "Score Breakdown — Top 20 Prioritized Targets",
         fontsize = 12, fontsize_number = 11)
dev.off()


write.csv(scored2, "Phase8_Target_Prioritization_Scores.csv", row.names = FALSE)

print(head(scored$symbol, 3))


# ── Druggability filter ───────────────────────────────────────────────────────
# Define stromal/ECM genes to penalize
stromal_genes <- c("MMP1", "MMP2", "MMP3", "MMP7", "MMP9", "MMP10", "MMP11", 
                   "MMP13", "MMP14", "COL1A1", "COL1A2", "COL3A1", "COL5A1", 
                   "COL5A2", "COL10A1", "COL11A1", "COL12A1", "FN1", "POSTN", 
                   "VCAN", "BGN", "AEBP1", "THBS2", "SPP1", "SULF1", "SULF2")

# Define known druggable/oncogenic genes to bonus
druggable_genes <- c("PARP1", "ERBB2", "KPNA2", "STMN1", "TUBA1C", "SHCBP1",
                     "AURKB", "PLK1", "CDK1", "BUB1", "CCNB1", "CCNA2",
                     "TOP2A", "TYMS", "RRM2", "PCNA", "MCM2", "MCM6",
                     "STAT1", "STAT3", "E2F1", "MYBL2", "FOXM1", "KIF11",
                     "KIF2C", "KIF20A", "CENPE", "MAD2L1", "BIRC5")

print(head(scored2$symbol, 3))

# ── Save updated scores ───────────────────────────────────────────────────────
write.csv(scored2, "Phase8_Target_Prioritization_Scores_v2.csv", row.names = FALSE)

