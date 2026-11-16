# ==============================================================================
# run_all.R - runs the full pipeline in order and records the session.
# Heavy steps (03, 04) use 1,000 CV repeats and 5,000 permutations: run them on a
# multi-core machine. Each script can also be run on its own.
# ==============================================================================
source(here::here("R", "00_setup.R"))

steps <- c("01_phenotypes.R", "02_preprocess_rnaseq.R", "03_splsda.R",
           "04_spls_carcass.R",
           "06_enrichment_overlap.R", "07_figure_stability.R")
for (s in steps) {
  message("\n>>> ", s)
  source(here::here("R", s), echo = FALSE)
}
writeLines(capture.output(sessionInfo()), file.path(dir_res, "sessionInfo.txt"))
