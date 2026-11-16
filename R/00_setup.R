# ==============================================================================
# 00_setup.R
# Global configuration, paths and helper functions. Sourced by every script.
# ==============================================================================

required <- c(
  "here", "readr", "dplyr", "tidyr", "tibble", "ggplot2",
  "lme4", "lmerTest", "nlme", "emmeans", "performance", "simr", "TOSTER",
  "edgeR", "limma", "DESeq2", "vsn",
  "mixOmics", "BiocParallel",
  "enrichR", "eulerr", "fgsea", "msigdbr"
)
missing_pkgs <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  stop("Missing packages: ", paste(missing_pkgs, collapse = ", "),
       "\nCRAN: install.packages(); Bioconductor: BiocManager::install().")
}

# ---- reproducibility ---------------------------------------------------------
SEED <- 1234
set.seed(SEED)

# ---- paths -------------------------------------------------------------------
dir_raw   <- here::here("data", "raw")        # count matrix (not tracked by git)
dir_meta <- here::here("data", "metadata")    # sample + phenotype tables
dir_proc <- here::here("data", "processed") # VST matrices etc.
dir_res   <- here::here("results")
dir_fig   <- here::here("figures")
for (d in c(dir_proc, dir_res, dir_fig)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---- analysis parameters (as described in Materials and Methods) --------------
N_CORES        <- max(1, parallel::detectCores() - 1)    # set to 1 on Windows
NFOLDS         <- 3        # M-fold cross-validation
NREPEAT_CV     <- 1000     # CV repeats (tuning, perf, stability)
NCOMP_MAX      <- 5        # components evaluated during tuning
N_PERM         <- 5000     # label permutations (as stated in the Methods)
# CV repeats INSIDE each permutation (3-fold CV, same structure as the observed model).
# 1,000 repeats x 5,000 permutations is computationally prohibitive, so each permuted dataset is
# evaluated with 10 repeats (stated in the Methods). Observed BER: NREPEAT_CV repeats.
NREPEAT_PERM   <- 10
STAB_THRESHOLD <- 0.50     # selection-frequency threshold
BOOT_B         <- 1000     # bootstrap replicates (sPLS)

# keepX grid for tune.splsda(): 1-100 in steps of 1 (the BER curves of Supplementary Figures
# S3c and S4c are evaluated at every integer up to 100). It contains the keepX of the reported
# final models (weaning 59/9; pre-slaughter 63/8). This is computationally heavy: reduce
# NREPEAT_CV for a quick test run.
KEEPX_GRID <- 1:100

# keepX of the final models reported in the manuscript (Results). tune.splsda() is still run;
# 03_splsda.R warns if its choice differs from these values.
KEEPX_PUBLISHED <- list(weaning = c(59, 9), preslaughter = c(63, 8))
stopifnot(all(unlist(KEEPX_PUBLISHED) %in% KEEPX_GRID))

# Number of genes after filtering reported in the manuscript (checked in script 02).
EXPECTED_GENES <- c(weaning = 12385, preslaughter = 11958)

# ---- input validation ---------------------------------------------------------
check_cols <- function(df, cols, file) {
  miss <- setdiff(cols, names(df))
  if (length(miss) > 0)
    stop(sprintf("%s is missing column(s): %s\nFound: %s", file,
                 paste(miss, collapse = ", "), paste(names(df), collapse = ", ")), call. = FALSE)
  invisible(df)
}

load_metadata <- function() {
  f <- file.path(dir_meta, "sample_metadata.csv")
  m <- readr::read_csv(f, show_col_types = FALSE)
  check_cols(m, c("sample_id", "animal_id", "treatment", "timepoint"), f)
  m$treatment <- factor(m$treatment, levels = c("CW", "EW"))
  m$timepoint <- factor(m$timepoint, levels = c("weaning", "preslaughter"))
  if (anyNA(m$treatment) || anyNA(m$timepoint))
    stop("sample_metadata.csv: treatment must be CW/EW and timepoint weaning/preslaughter.")
  m
}

load_counts <- function() {
  as.matrix(read.delim(file.path(dir_raw, "counts.tsv"), row.names = 1, check.names = FALSE))
}

load_phenotypes <- function() {
  f <- file.path(dir_meta, "phenotypes.csv")
  p <- readr::read_csv(f, show_col_types = FALSE)
  check_cols(p, c("animal_id", "treatment", "paddock", "rnaseq", "adg_120_205", "adg_205_914",



                    "adg_120_914", "bw_120", "bw_205", "bw_914", "hcw", "cy", "rea", "sft", "lipids"), f)
    p$treatment <- factor(p$treatment, levels = c("CW", "EW"))
    p
}

# ---- version-robust wrappers ---------------------------------------------------
# msigdbr changed its API in v10 (category/subcategory -> collection/subcollection).
get_go_bp_sets <- function(species = "Bos taurus") {
  d <- tryCatch(
    msigdbr::msigdbr(species = species, collection = "C5", subcollection = "GO:BP"),
    error = function(e) msigdbr::msigdbr(species = species, category = "C5", subcategory = "GO:BP"))
  split(d$gene_symbol, d$gs_name)
}

# Enrichr library names change between releases: resolve the newest library matching a
# pattern (e.g. "GO_Biological_Process") from the live list.
resolve_enrichr_dbs <- function(patterns = c("^GO_Biological_Process_20\\d\\d$",
                                             "^KEGG_20\\d\\d_Human$",
                                             "^Reactome_20\\d\\d$")) {
  avail <- enrichR::listEnrichrDbs()$libraryName
  vapply(patterns, function(p) {
    hit <- sort(grep(p, avail, value = TRUE), decreasing = TRUE)
    if (length(hit) == 0) stop("No Enrichr library matches ", p)
    hit[1]
  }, character(1), USE.NAMES = FALSE)
}
