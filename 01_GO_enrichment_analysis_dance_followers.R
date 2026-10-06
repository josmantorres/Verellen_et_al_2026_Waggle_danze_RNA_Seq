# ============================================================
# GO enrichment analysis using ermineR
# ============================================================
#
# Description:
# Gene Ontology (GO) enrichment analysis of RNA-seq differential
# expression results using the ermineR package and the Gene Score
# Resampling (GSR) method.
#
# The analysis is performed separately for:
#   - Biological Process (BP)
#   - Molecular Function (MF)
#   - Cellular Component (CC)
#
# Comparisons:
#   - Pollen vs Control
#   - Nectar vs Control
#   - Nectar vs Pollen
#   - Followers vs Non-followers
#
# Requirements:
#   - R
#   - ermineR
#   - rJava
#
# ============================================================


# ------------------------------------------------------------
# 1. Install required packages
# ------------------------------------------------------------

# Install ermineR from GitHub if necessary:
# install.packages("remotes")
# remotes::install_github("PavlidisLab/ermineR")

# rJava installation on Ubuntu may require system-level
# dependencies. See:
# https://www.r-bloggers.com/2018/02/installing-rjava-on-ubuntu/


# ------------------------------------------------------------
# 2. Load required packages
# ------------------------------------------------------------

library(ermineR)
library(rJava)
library(httr)


# Increase the timeout for downloading annotation resources
# when internet connection speed is slow.
options(timeout = max(10000, getOption("timeout")))


# ------------------------------------------------------------
# 3. Input and output directories
# ------------------------------------------------------------

# Replace these paths with the directories used on your system.
input_dir <- "path/to/input/files"
output_dir <- "path/to/output/files"


# ------------------------------------------------------------
# 4. Load GO annotation data
# ------------------------------------------------------------

# The GO annotation file must contain:
#   - a gene identifier column ("names")
#   - a GO term column ("Goterm")
#
# Each row represents a gene-GO term association.

goAnnots <- read.table(
  file.path(input_dir, "GO_terms_2024.txt"),
  sep = "\t",
  header = TRUE,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 5. Prepare GO annotation list for ermineR
# ------------------------------------------------------------

# Create a list of GO terms associated with each gene.
GOList <- split(
  goAnnots$Goterm,
  goAnnots$names
)

# Convert the GO annotation list to the format required by ermineR.
annotationList <- makeAnnotation(
  GOList,
  return = TRUE
)


# ------------------------------------------------------------
# 6. Load differential expression scores
# ------------------------------------------------------------

# The score files contain gene-level raw p-values from the
# differential expression analyses.
#
# Because ermineR is configured below with logTrans = TRUE and
# bigIsBetter = FALSE, smaller p-values receive higher scores
# after transformation.

scoreList_PvsC <- read.table(
  file.path(input_dir, "IDs_vs_pvalue_PvsC.txt"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

scoreList_NvsC <- read.table(
  file.path(input_dir, "IDs_vs_pvalue_NvsC.txt"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

scoreList_NvsP <- read.table(
  file.path(input_dir, "IDs_vs_pvalue_NvsP.txt"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

scoreList_FvsNF <- read.table(
  file.path(input_dir, "IDs_vs_pvalue_Followers_vs_Non-followers.txt"),
  sep = "\t",
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------
# 7. GO enrichment analysis
# ------------------------------------------------------------
#
# ermineR parameters:
#
#   aspects:
#       "B" = Biological Process
#       "M" = Molecular Function
#       "C" = Cellular Component
#
#   test = "GSR":
#       Gene Score Resampling
#
#   pAdjust = "FDR":
#       False Discovery Rate correction
#
#   iterations = 200000:
#       Number of resampling iterations
#
#   geneReplicates = "mean":
#       Mean score used for replicated gene identifiers
#
#   stats = "mean":
#       Mean statistic used by the GSR analysis
#
#   minClassSize = 10:
#       Minimum number of genes in a GO category
#
#   maxClassSize = 200:
#       Maximum number of genes in a GO category
#
#   logTrans = TRUE:
#       Log-transform the input scores
#
#   bigIsBetter = FALSE:
#       Smaller input p-values correspond to stronger evidence
#       for differential expression.
#
# ------------------------------------------------------------


# ============================================================
# Pollen vs Control
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "B",
  scores = scoreList_PvsC,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_BP_PvsC.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "M",
  scores = scoreList_PvsC,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_MF_PvsC.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "C",
  scores = scoreList_PvsC,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_CC_PvsC.txt"
  )
)


# ============================================================
# Nectar vs Control
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "B",
  scores = scoreList_NvsC,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_BP_NvsC.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "M",
  scores = scoreList_NvsC,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_MF_NvsC.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "C",
  scores = scoreList_NvsC,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_CC_NvsC.txt"
  )
)


# ============================================================
# Nectar vs Pollen
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "B",
  scores = scoreList_NvsP,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_BP_NvsP.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "M",
  scores = scoreList_NvsP,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_MF_NvsP.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "C",
  scores = scoreList_NvsP,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_CC_NvsP.txt"
  )
)


# ============================================================
# Followers vs Non-followers
# ============================================================

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "B",
  scores = scoreList_FvsNF,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_BP_FvsNF.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "M",
  scores = scoreList_FvsNF,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_MF_FvsNF.txt"
  )
)

ermineR(
  annotation = annotationList,
  expression = NULL,
  aspects = "C",
  scores = scoreList_FvsNF,
  scoreColumn = 1,
  logTrans = TRUE,
  bigIsBetter = FALSE,
  test = "GSR",
  pAdjust = "FDR",
  iterations = 200000,
  geneReplicates = "mean",
  stats = "mean",
  return = FALSE,
  minClassSize = 10,
  maxClassSize = 200,
  output = file.path(
    output_dir,
    "Enriched_GOterms_CC_FvsNF.txt"
  )
)
