# =========================================================
# 01. Differential expression analysis: honey bee dance followers
# =========================================================
# Experimental groups (column "treatment" of the metadata):
#   - Control : non-followers
#   - Nectar  : followers of nectar foragers
#   - Pollen  : followers of pollen foragers
#
# Comparisons (all with the SAME pipeline):
#   1) Followers (Nectar + Pollen) vs Control   [main comparison]
#   2) Nectar vs Control
#   3) Pollen vs Control
#   4) Pollen vs Nectar
#
# Model:      negative binomial GLM, Wald test (H0: log2FC = 0)
# Shrinkage:  apeglm (lfcShrink)
# DEG rule:   apeglm s-value < 0.05 (no log2FC magnitude threshold)
# Settings:   cooksCutoff = TRUE; independentFiltering = FALSE
#             (independent filtering only affects adjusted p-values and
#             does not influence s-value-based calls)
#
# No gene-filtering step is applied before the analysis; genes with
# zero counts across all samples are excluded automatically by DESeq2.
#
#
# Usage (from the repository root):
#   Rscript scripts/01_DE_dance_followers.R
#
# Input files (tab-delimited, in `input/`; see README for details):
#   counts_matrix.txt   : gene IDs in column 1, raw integer counts after
#   sample_metadata.txt : columns "sample" and "treatment"
#                         (Control / Nectar / Pollen)
# =========================================================

library(DESeq2)
library(apeglm)
library(pheatmap)
library(RColorBrewer)

# ---------------------------------------------------------
# 0. Parameters
# ---------------------------------------------------------

input_dir  <- "input"
output_dir <- "output"
alpha_s    <- 0.05    # s-value threshold (also alpha for results/summary)
independent_filtering <- FALSE
dec_out    <- "."     # decimal separator of exported tables

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
out <- function(x) file.path(output_dir, x)

# ---------------------------------------------------------
# 1. Import counts and metadata
# ---------------------------------------------------------

counts <- read.delim(file.path(input_dir, "counts_matrix.txt"),
                     sep = "\t", header = TRUE)
count_matrix <- data.matrix(counts[, 2:ncol(counts)])
rownames(count_matrix) <- counts[, 1]

meta <- read.delim(file.path(input_dir, "sample_metadata.txt"),
                   sep = "\t", header = TRUE)
stopifnot(all(c("sample", "treatment") %in% colnames(meta)))

# Samples in metadata must match the count columns (same order)
if (!all(colnames(count_matrix) == meta$sample)) {
  stop("Sample names/order differ between the counts and the metadata.")
}
rownames(meta) <- meta$sample

stopifnot(all(meta$treatment %in% c("Control", "Nectar", "Pollen")))
meta$treatment <- factor(meta$treatment,
                         levels = c("Control", "Nectar", "Pollen"))

# Followers (Nectar + Pollen) vs Non-followers (Control)
meta$follow_state <- factor(
  ifelse(meta$treatment == "Control", "NonFollower", "Follower"),
  levels = c("NonFollower", "Follower")
)

# ---------------------------------------------------------
# 2. Exploratory analyses (all genes, blind rlog)
# ---------------------------------------------------------

dds_explore <- DESeqDataSetFromMatrix(count_matrix, meta,
                                      design = ~ treatment)
rld <- rlogTransformation(dds_explore, blind = TRUE, fitType = "parametric")

write.table(data.frame(assay(rld)), out("rlog_values.txt"),
            sep = "\t", dec = dec_out,
            row.names = TRUE, col.names = NA, quote = FALSE)

sample_dists <- dist(t(assay(rld)))                # Euclidean distance
dist_matrix  <- as.matrix(sample_dists)

tiff(out("Sample_distance_heatmap.tif"), width = 12, height = 12,
     units = "cm", res = 300, pointsize = 8, compression = "lzw")
pheatmap(dist_matrix,
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
dev.off()

# PCA on the 500 most variable genes (plotPCA default: ntop = 500)
tiff(out("PCA_rlog.tif"), width = 15, height = 15, units = "cm",
     res = 300, pointsize = 10, compression = "lzw")
print(plotPCA(rld, intgroup = "treatment"))
dev.off()

# ---------------------------------------------------------
# 3. Common function applied to ALL comparisons
# ---------------------------------------------------------
# - results(): Wald test of H0: LFC = 0; cooksCutoff = TRUE sets
#   pvalue = NA for genes flagged as outliers by Cook's distance.
# - lfcShrink(type = "apeglm"): requires a model coefficient (not a
#   contrast); svalue = TRUE is needed to obtain s-values. The s-value is
#   the estimated probability that the sign of the shrunken log2 fold
#   change is incorrect (false sign rate).
# - apeglm does not apply Cook's cutoff by itself: the s-value is set to
#   NA for genes flagged as outliers by results().

run_apeglm <- function(dds, coef) {
  stopifnot(coef %in% resultsNames(dds))

  res <- results(dds,
                 name = coef,
                 alpha = alpha_s,
                 independentFiltering = independent_filtering,
                 cooksCutoff = TRUE)

  res_shr <- lfcShrink(dds,
                       coef = coef,
                       res = res,
                       type = "apeglm",
                       svalue = TRUE)

  res_shr$svalue[is.na(res$pvalue)] <- NA

  # Wald statistics (unshrunken estimates, H0: LFC = 0)
  res_shr$pvalue_Wald <- res$pvalue
  res_shr$padj_Wald   <- res$padj

  res_shr
}

# ---------------------------------------------------------
# 4. Model 1: Followers vs Non-followers
# ---------------------------------------------------------

dds_follow <- DESeqDataSetFromMatrix(count_matrix, meta,
                                     design = ~ follow_state)
dds_follow$follow_state <- relevel(dds_follow$follow_state, ref = "NonFollower")
dds_follow <- DESeq(dds_follow, test = "Wald")

resultsNames(dds_follow)
# "Intercept" "follow_state_Follower_vs_NonFollower"

res_FvsC <- run_apeglm(dds_follow, "follow_state_Follower_vs_NonFollower")

# ---------------------------------------------------------
# 5. Model 2: resource-specific comparisons (reference = Control)
# ---------------------------------------------------------

dds_treatment <- DESeqDataSetFromMatrix(count_matrix, meta,
                                        design = ~ treatment)
dds_treatment$treatment <- relevel(dds_treatment$treatment, ref = "Control")
dds_treatment <- DESeq(dds_treatment, test = "Wald")

resultsNames(dds_treatment)
# "Intercept" "treatment_Nectar_vs_Control" "treatment_Pollen_vs_Control"

res_NvsC <- run_apeglm(dds_treatment, "treatment_Nectar_vs_Control")
res_PvsC <- run_apeglm(dds_treatment, "treatment_Pollen_vs_Control")

# ---------------------------------------------------------
# 6. Model 3: Pollen vs Nectar (reference = Nectar)
# ---------------------------------------------------------
# apeglm does not accept contrasts: the model is refitted with Nectar as
# the reference level so that Pollen vs Nectar is a model coefficient.

dds_PN <- dds_treatment
dds_PN$treatment <- relevel(dds_PN$treatment, ref = "Nectar")
dds_PN <- DESeq(dds_PN, test = "Wald")

resultsNames(dds_PN)
# "Intercept" "treatment_Control_vs_Nectar" "treatment_Pollen_vs_Nectar"

res_PvsN <- run_apeglm(dds_PN, "treatment_Pollen_vs_Nectar")

