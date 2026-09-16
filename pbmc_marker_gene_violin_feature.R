#!/usr/bin/env Rscript
# =============================================================================
# PBMC marker-gene violin + feature plots (B, pDC, T_CD4, T_CD8, Mono only)
# -----------------------------------------------------------------------------
# RNA expression (from the aggregated Cut&Tag+DynaTag combined object) for a
# curated marker-gene panel, restricted to 5 major groups: B, pDC, T_CD4,
# T_CD8, Mono. One INDIVIDUAL violin plot and one INDIVIDUAL feature plot per
# gene (not faceted grids), styled for presentation: big titles/axis/legend
# text, same violin+boxplot+median-label convention used throughout this
# session's QC violins (theme_classic), feature plots on the existing
# "umap_naive" RNA embedding (subset to these 5 clusters), theme_minimal.
#
# Gene panel (protein/CD markers translated to gene symbols; user-specified,
# 2026-09-04). pDC had no genes specified -- CLEC4C/IRF7 added as canonical
# pDC markers since pDC is one of the 5 requested clusters.
#   B:      CD19, MS4A1 (CD20), IGHM (IGH), CD27 (memory), CD38 (plasma)
#   T:      CD3D, CD3E (CD3), IL2
#   T_CD4:  CD4, SELL (CD62L), IL7R (CD127), HLA-DRA, HLA-DRB1 (HLA-DR), IFNG, TNF
#   T_CD8:  CD8A, CD8B, PRF1 (Perforin), GZMB (Granzyme), TRAC, TRBC1, TRBC2 (TCR a/b)
#   Mono:   CD14, FCGR3A (CD16), ITGAM (CD11b), CD33
#   DC:     FCGR3A (CD16), ITGAM (CD11b), ITGAX (CD11c)
#   pDC:    CLEC4C, IRF7 (added -- none given)
#
# PLUS the top-5 FindAllMarkers hits per cluster from pbmc_marker_gene_heatmap.R
# (PBMCs_marker_genes_top.csv), for these same 5 clusters, de-duplicated
# against the curated list above.
#
# Run:
#   source /home/mattia/miniconda3_n/etc/profile.d/conda.sh
#   conda run -n nanoctarna-analysis Rscript scripts/pbmc_marker_gene_violin_feature.R
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2)
})

obj_dir <- "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/first_eval/objects"
fig_dir <- "/date/gcb/gcb_MZ/multiNanoCT/samples/PBMCs-MS/figures_manual/marker_genes"
dir.create(file.path(fig_dir, "violin"),  recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(fig_dir, "feature"), recursive = TRUE, showWarnings = FALSE)

CLUSTERS <- c("B", "pDC", "T_CD4", "T_CD8", "Mono")
CLUSTER_COLORS <- c(B = "#1f77b4", pDC = "#d62728", T_CD4 = "#2ca02c",
                    T_CD8 = "#17becf", Mono = "#ff7f0e")

GENES_CURATED <- c("GRCh38-CXCR3","GRCh38-TBX21", "GRCh38-MS4A1")                                    # pDC (added)

# top-5 FindAllMarkers hits (pbmc_marker_gene_heatmap.R) for these same 5
# clusters (PBMCs_marker_genes_top.csv), in cluster order; de-duplicated
# against the curated list above (MS4A1/CD8A/IL7R already present).
GENES_DISCOVERED <- c(
  "FCRL1", "EBF1", "LINC00926", "BANK1",              # B (MS4A1 already curated)
  "MSR1", "MS4A7", "CD300E", "CLEC7A", "AQP9",        # Mono
  "TSHZ2", "FAAH2", "ICOS", "SESN3", "INPP4B",        # T_CD4
  "SERPINF1", "RHEX", "PALD1", "VASH2", "LINC01478",  # pDC
  "SGCD", "TRGC2", "CCL5", "CADM1"                    # T_CD8 (CD8A already curated)
)

GENES <- unique(c(GENES_CURATED, GENES_DISCOVERED))

cat("Loading combined (aggregated) object...\n")
combined <- readRDS(file.path(obj_dir, "PBMCs_peaks_combined_naive_harmony.rds"))
rna <- combined[["PBMC|RNA"]]
DefaultAssay(rna) <- "RNA"

missing <- setdiff(GENES, rownames(rna))
if (length(missing)) cat("WARNING missing genes (skipped):", paste(missing, collapse=", "), "\n")
GENES <- intersect(GENES, rownames(rna))

sub <- subset(rna, subset = predicted.id_scRNA_Atlas %in% CLUSTERS)
sub$predicted.id_scRNA_Atlas <- factor(sub$predicted.id_scRNA_Atlas, levels = CLUSTERS)
Idents(sub) <- "predicted.id_scRNA_Atlas"
cat(sprintf("Subset to %d cells across %s\n", ncol(sub), paste(CLUSTERS, collapse=", ")))
print(table(sub$predicted.id_scRNA_Atlas))

emb <- Embeddings(sub, "umap_naive")

## ---- individual violin plot (theme_classic, big fonts, median label) ------
gene_violin <- function(gene) {
  d <- data.frame(expr = FetchData(sub, vars = gene)[, 1],
                  cluster = sub$predicted.id_scRNA_Atlas)
  ggplot(d, aes(x = cluster, y = expr, fill = cluster)) +
    geom_violin(scale = "width", trim = TRUE, alpha = 0.85) +
    geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.85) +
    stat_summary(fun = median, geom = "text",
                 aes(label = round(after_stat(y), 2)), vjust = -0.6, size = 6) +
    scale_fill_manual(values = CLUSTER_COLORS) +
    ylab(paste0(gene, " expression (log-normalized)")) + xlab("") +
    ggtitle(gene) +
    theme_classic(base_size = 20) +
    theme(legend.position = "none",
          plot.title  = element_text(size = 28, face = "bold.italic", hjust = 0.5),
          axis.title.y = element_text(size = 20),
          axis.text.x  = element_text(size = 22, angle = 45, hjust = 1),
          axis.text.y  = element_text(size = 16))
}

## ---- individual feature plot (theme_minimal, big fonts) --------------------
gene_feature <- function(gene) {
  d <- data.frame(UMAP_1 = emb[,1], UMAP_2 = emb[,2],
                  expr = FetchData(sub, vars = gene)[, 1])
  d <- d[order(d$expr), ]  # plot highest-expressing cells last (on top)
  ggplot(d, aes(UMAP_1, UMAP_2, color = expr)) +
    geom_point(size = 1.1, alpha = 0.9) +
    scale_color_gradient(low = "lightgrey", high = "#0033CC", name = "Expr.") +
    ggtitle(gene) +
    theme_minimal(base_size = 20) +
    theme(plot.title = element_text(size = 28, face = "bold.italic", hjust = 0.5),
          axis.title = element_text(size = 20),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          legend.text  = element_text(size = 16),
          legend.title = element_text(size = 18),
          panel.grid = element_blank())
}

save_both <- function(p, path, w, h) {
  for (ext in c("png", "pdf")) ggsave(paste0(path, ".", ext), p, width = w, height = h, dpi = 300)
}

for (gene in GENES) {
  cat("  plotting:", gene, "\n")
  save_both(gene_violin(gene),  file.path(fig_dir, "violin",  paste0("Vln_",  gsub("[^A-Za-z0-9]", "_", gene))), w = 6.5, h = 6.5)
  save_both(gene_feature(gene), file.path(fig_dir, "feature", paste0("Feat_", gsub("[^A-Za-z0-9]", "_", gene))), w = 7,   h = 6.5)
}

cat(sprintf("\nWrote %d violin + %d feature plots to %s\n", length(GENES), length(GENES), fig_dir))
cat("Done.\n")
