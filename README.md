# Differential expression in honey bee dance followers

This repository contains the code used for trimming and mapping/counting reads, the read counts matrix, the experimental design information, the R code used for the differential expression and GO-enrichment analyses and the output files, Rlog-transformed values, and TPM tables in the manuscript:

> [Full citation of the manuscript: authors, title, journal, year, DOI]

[![DOI](https://zenodo.org/badge/DOI/[ZENODO-DOI].svg)](https://doi.org/[ZENODO-DOI])

## Overview

The script analyzes RNA-seq data from honey bees belonging to three experimental groups:

- **Control:** non-followers of the waggle dance.
- **Nectar:** [followers of the waggle dance of nectar foragers].
- **Pollen:** [followers of the waggle dance of pollen foragers].

Gene-level raw counts are analyzed with DESeq2. **No gene-filtering step is applied**: genes with zero counts across all samples are excluded automatically by DESeq2.

### Comparisons

All four comparisons are run with the same pipeline:

1. **Followers (Nectar + Pollen) vs Control** (main comparison)
2. **Nectar vs Control**
3. **Pollen vs Control**
4. **Pollen vs Nectar**

Positive log2 fold changes indicate higher expression in the first group of each comparison.

### Statistical analysis

- **Model:** negative binomial GLM, Wald test (null hypothesis: log2FC = 0). No log2 fold-change magnitude threshold is applied.
  - Comparison 1 uses the model `~ follow_state` (Follower vs NonFollower).
  - Comparisons 2 and 3 use `~ treatment` with Control as the reference level.
  - Comparison 4 uses the same model with Nectar as the reference level (apeglm requires a model coefficient, so the model is refitted).
- **Effect-size shrinkage:** apeglm (`lfcShrink`).
- **DEG criterion:** apeglm s-value < 0.05. The s-value estimates the probability that the sign of the shrunken log2 fold change is incorrect (false sign rate). It is not a false discovery rate, and it does not impose a minimum magnitude of change.
- **Other settings:** `cooksCutoff = TRUE`. Because apeglm does not apply Cook's cutoff by itself, the s-value is set to `NA` for genes flagged as outliers. `independentFiltering = FALSE`; independent filtering only affects adjusted p-values, so it does not change s-value-based calls (it can be switched with the `independent_filtering` parameter).
- **Tables:** include the shrunken log2 fold change and s-value, plus the Wald p-value and adjusted p-value (both computed with the Wald test on unshrunken estimates, null hypothesis log2FC = 0; reported for information only).
- **Exploratory analyses:** rlog transformation (`blind = TRUE`), sample-distance heatmap (Euclidean distance) and PCA (500 most variable genes).


## Repository structure

```
.
├── scripts/
│   └── Trimming and Mapping/Counting code   # Code used for trimming the raw reads with Trimmomatic and Mapping/Count with STAR
│   └── 01_DE_dance_followers.R   # Code used for differential expression analysis with DESeq2
│   └── 01_GO_enrichment_analysis_dance_followers.R   # Code used for GO enrichment analysis with ermineR
├── input/                     # input files (see below)
│   └── Raw reads count matrix.txt
│   └── Experimental_design.txt   # Samples and treatments information
│   └── Raw_reads_info.txt   # BioProject code  
├── output/                    # created when the script is run
│   └── Deseq2_statistics_output_Followers_vs_Control.txt   # DESeq2 analysis output from the Followers vs Non-followers (Control) comparison
│   └── Deseq2_statistics_output_Nectar_vs_Control.txt   # DESeq2 analysis output from the Nectar vs Control (Non-followers) comparison
│   └── Deseq2_statistics_output_Pollen_vs_Control.txt   # DESeq2 analysis output from the Pollen vs Control (Non-followers) comparison
│   └── Deseq2_statistics_output_Pollen_vs_Nectar.txt   # DESeq2 analysis output from the Pollen vs Nectar comparison
│   └── GO_enrichment_BP_FvsNF.txt   # GO-term enrichment analysis output for biological process for Follower vs Control (Non-followers)
│   └── GO_enrichment_BP_NvsC.txt   # GO-term enrichment analysis output for biological process for Nectar vs Control (Non-followers)
│   └── GO_enrichment__BP_PvsC.txt   # GO-term enrichment analysis output for biological process for Pollen vs Control (Non-followers)
│   └── GO_enrichment__BP_PvsN.txt   # GO-term enrichment analysis output for biological process for Pollen vs Nectar
│   └── GO_enrichment__MF_FvsNF.txt   # GO-term enrichment analysis output for molecular function for Follower vs Control (Non-followers)
│   └── GO_enrichment_MF_NvsC.txt   # GO-term enrichment analysis output for molecular function for Nectar vs Control (Non-followers)
│   └── GO_enrichment__MF_PvsC.txt   # GO-term enrichment analysis output for molecular function for Pollen vs Control (Non-followers)
│   └── GO_enrichment__MF_PvsN.txt   # GO-term enrichment analysis output for molecular function for Pollen vs Nectar
│   └── GO_enrichment__CC_FvsNF.txt   # GO-term enrichment analysis output for cellular component for Follower vs Control (Non-followers)
│   └── GO_enrichment__CC_NvsC.txt   # GO-term enrichment analysis output for cellular component for Nectar vs Control (Non-followers)
│   └── GO_enrichment__CC_PvsC.txt   # GO-term enrichment analysis output for cellular component for Pollen vs Control (Non-followers)
│   └── GO_enrichment__CC_PvsN.txt   # GO-term enrichment analysis output for cellular component for Pollen vs Nectar
│   └── Experimental_design.txt   # Samples and treatments information
│   └── RLOG_Transformed_values.txt
│   └── TPM_values.txt   # Transcripts per million values table
└── README.md
```

## Requirements

- R [version 2026.09.0]
- R packages: `DESeq2` [1.42.1], `apeglm` [1.24], `pheatmap` [1.0.13], `RColorBrewer` [1-1.3]

The exact versions used are written to `output/sessionInfo.txt`.

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm"))
install.packages(c("pheatmap", "RColorBrewer"))
```

## Input files

Both files are tab-delimited with a header row and must be placed in `input/`.

**`input/Raw reads count matrix.txt`**: first column = gene IDs; remaining columns = raw (unnormalized) integer counts, one column per sample.

| gene_id | sample1 | sample2 | ... |
|---------|---------|---------|-----|
| gene_A  | 120     | 98      | ... |

**`input/Experimental_design.txt`**: one row per sample, with the same sample names and order as the count columns. Columns:

- `sample`: sample name.
- `treatment`: `Control`, `Nectar` or `Pollen`.

The script stops with an error message if the sample names or order do not match between the counts and the metadata.

## Usage

From the repository root:

```bash
Rscript scripts/01_DE_dance_followers.R
```

or open the script in RStudio (with the repository folder as working directory). Paths, thresholds and the decimal separator are set in the "Parameters" section at the top of the script.

## Outputs

| File | Description |
|------|-------------|
| `GO_enrichment_BP_FvsNF.txt` | GO-term enrichment analysis output for biological process for Follower vs Control (Non-followers)
| `GO_enrichment_BP_NvsC.txt` | GO-term enrichment analysis output for biological process for Nectar vs Control (Non-followers)
| `GO_enrichment_BP_PvsC.txt` | GO-term enrichment analysis output for biological process for Pollen vs Control (Non-followers) |
| `GO_enrichment_BP_PvsN.txt` | GO-term enrichment analysis output for biological process for Pollen vs Nectar |
| `GO_enrichment_MF_FvsNF.txt` | GO-term enrichment analysis output for molecular function for Follower vs Control (Non-followers) |
| `GO_enrichment_MF_NvsC.txt` | GO-term enrichment analysis output for molecular function for Nectar vs Control (Non-followers) |
| `GO_enrichment_MF_PvsC.txt` | GO-term enrichment analysis output for molecular function for Pollen vs Control (Non-followers) |
| `GO_enrichment_MF_PvsN.txt` | GO-term enrichment analysis output for molecular function for Pollen vs Nectar |
| `GO_enrichment_CC_FvsNF.txt` | GO-term enrichment analysis output for cellular component for Follower vs Control (Non-followers) |
| `GO_enrichment_CC_NvsC.txt` | GO-term enrichment analysis output for cellular component for Nectar vs Control (Non-followers) |
| `GO_enrichment_CC_PvsC.txt` | GO-term enrichment analysis output for cellular component for Pollen vs Control (Non-followers) |
| `GO_enrichment_CC_PvsN.txt` | GO-term enrichment analysis output for cellular component for Pollen vs Nectar |
| `Experimental_design.txt` | Samples and treatments information |
| `RLOG_Transformed_values.txt` | 
| `TPM_values.txt` | Transcripts per million values table |

`<comparison>` is one of `Followers_vs_Control`, `Nectar_vs_Control`, `Pollen_vs_Control` and `Pollen_vs_Nectar`.

## Data availability

- Raw sequencing data: NCBI SRA: SRR34735861; SRR34735864; SRR34735863; SRR34735862; SRR34735860; SRR34735865
- Count table used in this analysis: included in `input/Raw reads count matrix.txt` 

## Citation

If you use this code, please cite the manuscript above and the archived version of this repository:

> [Authors]. [Repository title]. Zenodo. [Year]. https://doi.org/[ZENODO-DOI]

DESeq2: Love MI, Huber W, Anders S (2014). Genome Biology 15:550.
apeglm: Zhu A, Ibrahim JG, Love MI (2019). Bioinformatics 35:2084–2092.

## License

[e.g., MIT License; see `LICENSE`]

## Contact

[Name, institution, e-mail]
