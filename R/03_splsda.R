# ==============================================================================
# 03_splsda.R
# sPLS-DA (EW vs CW) at each time point with mixOmics:
#    tuning (ncomp, keepX) by repeated M-fold CV -> final model -> perf() ->
#    selection-frequency (stability) -> label-permutation test.
# NOTE: mixOmics has no built-in permutation test for sPLS-DA, so the permutation
#        is implemented explicitly below (shuffle Y, rerun the same pipeline, record BER).
# Inputs : data/processed/vst_<timepoint>.rds, data/metadata/sample_metadata.csv
# Outputs: results/splsda_<tp>_summary.csv
#           results/splsda_<tp>_stability.csv      (frequency per gene and component)
#           results/splsda_<tp>_selected.csv       (final-model genes + loadings)
#           results/splsda_<tp>_permutation.rds
#           results/splsda_<tp>_model.rds, figures/splsda_<tp>_*.pdf
# ==============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
   library(mixOmics); library(BiocParallel); library(dplyr); library(readr); library(tibble)
})

# ---- permutation test ----------------------------------------------------------
permutation_ber <- function(X, Y, ncomp, keepX, n_perm = N_PERM,
                            nrepeat = NREPEAT_PERM, folds = NFOLDS) {
  old_kind <- RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG"); set.seed(SEED + 1)
  on.exit(RNGkind(old_kind), add = TRUE)
  one <- function(i) {
    Yp <- sample(Y)
    m <- mixOmics::splsda(X, Yp, ncomp = ncomp, keepX = keepX)
    p <- mixOmics::perf(m, validation = "Mfold", folds = folds,
                         nrepeat = nrepeat, progressBar = FALSE)
    p$error.rate$BER[ncomp, "max.dist"]
  }
  unlist(parallel::mclapply(seq_len(n_perm), one, mc.cores = N_CORES, mc.set.seed = TRUE))
}

run_splsda <- function(tp, keepx_published = NULL) {
  message("== sPLS-DA: ", tp)
  vst <- readRDS(file.path(dir_proc, sprintf("vst_%s.rds", tp)))
  meta <- load_metadata() |> dplyr::filter(timepoint == tp)
  meta <- meta[match(colnames(vst), meta$sample_id), ]
  X <- t(vst)
  Y <- droplevels(meta$treatment)
  bp <- BiocParallel::MulticoreParam(workers = N_CORES, RNGseed = SEED)

  # 1) tuning: number of components and keepX (BER, max.dist)
  tune <- mixOmics::tune.splsda(X, Y, ncomp = NCOMP_MAX, test.keepX = KEEPX_GRID,
                                validation = "Mfold", folds = NFOLDS,
                                nrepeat = NREPEAT_CV, dist = "max.dist",
                                measure = "BER", BPPARAM = bp, progressBar = FALSE)
  pdf(file.path(dir_fig, sprintf("splsda_%s_tuning.pdf", tp)), width = 6, height = 4.5)
  print(plot(tune))
  dev.off()

  ncomp <- tune$choice.ncomp$ncomp
  keepX <- tune$choice.keepX[seq_len(ncomp)]
  if (!is.null(keepx_published)) {
    if (!identical(as.numeric(keepX), as.numeric(keepx_published)))
      warning("keepX from tuning (", paste(keepX, collapse = ","), ") differs from the published ",
              "values (", paste(keepx_published, collapse = ","), "); using the published values.")
    keepX <- keepx_published
    ncomp <- length(keepX)
  }

  # 2) final model and cross-validated performance
  final <- mixOmics::splsda(X, Y, ncomp = ncomp, keepX = keepX)
  perf_res <- mixOmics::perf(final, validation = "Mfold", folds = NFOLDS,
                             nrepeat = NREPEAT_CV, auc = TRUE,
                             BPPARAM = bp, progressBar = FALSE)
  ber_obs <- perf_res$error.rate$BER[ncomp, "max.dist"]

  # 3) selection frequency (stability) per component
  stab <- bind_rows(lapply(seq_len(ncomp), function(h) {
    s <- perf_res$features$stable[[h]]
    tibble(gene = names(s), frequency = as.numeric(s), component = paste0("Comp", h))
  })) |> arrange(component, desc(frequency))
  write_csv(stab, file.path(dir_res, sprintf("splsda_%s_stability.csv", tp)))

  # 4) genes in the final model with loadings
  sel <- bind_rows(lapply(seq_len(ncomp), function(h) {
    v <- mixOmics::selectVar(final, comp = h)$value
    tibble(gene = rownames(v), loading = v$value.var, component = paste0("Comp", h))
  })) |>
    left_join(stab, by = c("gene", "component")) |>
    mutate(robust = frequency >= STAB_THRESHOLD,
           # direction from the data (the sign of mixOmics loadings is arbitrary)
           mean_EW = rowMeans(vst[gene, Y == "EW", drop = FALSE]),
           mean_CW = rowMeans(vst[gene, Y == "CW", drop = FALSE]),



             diff_EW_minus_CW = mean_EW - mean_CW,
             direction = ifelse(diff_EW_minus_CW > 0, "Higher in EW", "Higher in CW"))
    write_csv(sel, file.path(dir_res, sprintf("splsda_%s_selected.csv", tp)))

    # 5) permutation test
    perm_ber <- permutation_ber(X, Y, ncomp, keepX)
    p_perm <- (sum(perm_ber <= ber_obs) + 1) / (length(perm_ber) + 1)
    saveRDS(list(observed = ber_obs, permuted = perm_ber, p = p_perm),
            file.path(dir_res, sprintf("splsda_%s_permutation.rds", tp)))
    pdf(file.path(dir_fig, sprintf("splsda_%s_permutation.pdf", tp)), width = 5, height = 4)
    hist(perm_ber, breaks = 30, col = "grey75", border = "white",
         main = sprintf("Permutation test (p = %.4g)", p_perm), xlab = "Balanced error rate (BER)")
    abline(v = ber_obs, col = "red", lwd = 2)
    dev.off()

    # 6) sample plot
    pdf(file.path(dir_fig, sprintf("splsda_%s_indiv.pdf", tp)), width = 5, height = 5)
    mixOmics::plotIndiv(final, comp = 1:2, group = Y, ellipse = TRUE, legend = TRUE,
                        title = sprintf("sPLS-DA - %s", tp))
    dev.off()

    summ <- tibble(
      timepoint = tp, ncomp = ncomp,
      keepX = paste(keepX, collapse = ";"),
      BER_final = ber_obs, BER_permuted_mean = mean(perm_ber), p_permutation = p_perm,
      n_perm = length(perm_ber), nrepeat_perm = NREPEAT_PERM,
      var_expl_comp = paste(round(100 * final$prop_expl_var$X[seq_len(ncomp)]), collapse = ";"),
      n_robust_comp1 = sum(sel$robust[sel$component == "Comp1"], na.rm = TRUE)
    )
    write_csv(summ, file.path(dir_res, sprintf("splsda_%s_summary.csv", tp)))
    saveRDS(final, file.path(dir_res, sprintf("splsda_%s_model.rds", tp)))
    summ
}

res <- bind_rows(
  run_splsda("weaning",      KEEPX_PUBLISHED$weaning),
  run_splsda("preslaughter", KEEPX_PUBLISHED$preslaughter)
)
print(res)
# Reference values in the manuscript figures: weaning BER = 0.103 (random 0.482, p = 0.01);
# pre-slaughter BER = 0.418 (random 0.468, p = 0.277). Compare with splsda_*_summary.csv.
