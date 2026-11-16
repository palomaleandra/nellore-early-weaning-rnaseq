# ==============================================================================
# 02_preprocess_rnaseq.R
# Low-count filtering (edgeR) and variance-stabilizing transformation (DESeq2),
# separately for each time point.
# Inputs : data/raw/counts.tsv (genes x samples, featureCounts, ARS-UCD1.2)
#           data/metadata/sample_metadata.csv
#           columns: sample_id, animal_id, treatment (CW/EW), timepoint (weaning/preslaughter)
# Outputs: data/processed/vst_<timepoint>.rds / .csv
#           results/filtering_summary.csv
#           figures/meanSd_<timepoint>.pdf, figures/vst_boxplot_<timepoint>.pdf
# ==============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
   library(edgeR); library(DESeq2); library(vsn); library(dplyr); library(readr)
})

counts <- load_counts()
meta   <- load_metadata()
stopifnot(all(meta$sample_id %in% colnames(counts)))

summary_rows <- list()
for (tp in levels(meta$timepoint)) {
  m <- meta |> dplyr::filter(timepoint == tp)
  cn <- counts[, m$sample_id, drop = FALSE]

  # keep genes with >= 10 reads in at least n/2 libraries of this time point
  keep <- rowSums(cn >= 10) >= ceiling(ncol(cn) / 2)
  cn_f <- cn[keep, , drop = FALSE]

  dds <- DESeqDataSetFromMatrix(countData = round(cn_f),
                                colData = data.frame(m, row.names = m$sample_id),
                                design = ~ treatment)
  vsd <- DESeq2::vst(dds, blind = TRUE)        # blind = TRUE: QC only, design-agnostic
  vst_mat <- SummarizedExperiment::assay(vsd)

  saveRDS(vst_mat, file.path(dir_proc, sprintf("vst_%s.rds", tp)))
  write.csv(vst_mat, file.path(dir_proc, sprintf("vst_%s.csv", tp)))

  # diagnostics (Supplementary Figures S3 and S4)
  pdf(file.path(dir_fig, sprintf("meanSd_%s.pdf", tp)), width = 5, height = 4)
  print(vsn::meanSdPlot(vst_mat))
  dev.off()
  pdf(file.path(dir_fig, sprintf("vst_boxplot_%s.pdf", tp)), width = 6, height = 4)
  boxplot(vst_mat, las = 2, cex.axis = 0.6, ylab = "VST counts")
  dev.off()

  summary_rows[[tp]] <- tibble(timepoint = tp, n_samples = ncol(cn),
                               genes_before = nrow(cn), genes_after = nrow(cn_f))
}
filt <- bind_rows(summary_rows)
write_csv(filt, file.path(dir_res, "filtering_summary.csv"))
filt$expected_in_text <- EXPECTED_GENES[filt$timepoint]
filt$matches_text <- filt$genes_after == filt$expected_in_text
write_csv(filt, file.path(dir_res, "filtering_summary.csv"))
print(filt)
if (any(!filt$matches_text))
  warning("Gene counts after filtering differ from the manuscript (see filtering_summary.csv): ",
          "correct the Results/Methods text (use the values from this table) or the filter.")
