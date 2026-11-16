# ==============================================================================
# 04_spls_carcass.R
# Integrative sPLS: pre-slaughter transcripts (X) vs. five standardized carcass
# traits (Y). Because cross-validated Q2 can be <= 0 with n = 16, the model is treated
# as a DESCRIPTIVE covariance analysis and complemented with:
#    - repeated-CV Q2 and R2 (reported as they are)
#    - bootstrap selection frequency and CIs of gene-trait correlations
#    - leave-one-animal-out retention of the selected genes
#    - permutation test of the first-component X/Y association
# Inputs : data/processed/vst_preslaughter.rds, data/metadata/{sample_metadata,phenotypes}.csv
# Outputs: results/spls_*.csv, results/spls_model.rds, figures/spls_*.pdf
# ==============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
   library(mixOmics); library(BiocParallel); library(dplyr); library(readr); library(tibble)
})

vst   <- readRDS(file.path(dir_proc, "vst_preslaughter.rds"))
meta <- load_metadata() |> dplyr::filter(timepoint == "preslaughter")
pheno <- load_phenotypes()
meta <- meta[match(colnames(vst), meta$sample_id), ]
pheno <- pheno[match(meta$animal_id, pheno$animal_id), ]

traits <- c("hcw", "cy", "rea", "sft", "lipids")
X <- t(vst)
Y <- scale(as.matrix(pheno[, traits]))            # z-score standardization
rownames(Y) <- rownames(X)
bp <- BiocParallel::MulticoreParam(workers = N_CORES, RNGseed = SEED)

# 1) number of components by Q2 --------------------------------------------------
pls0 <- mixOmics::pls(X, Y, ncomp = NCOMP_MAX, mode = "regression")
perf0 <- mixOmics::perf(pls0, validation = "Mfold", folds = NFOLDS, nrepeat = NREPEAT_CV,
                        BPPARAM = bp, progressBar = FALSE)
q2 <- perf0$measures$Q2.total$summary
write_csv(as.data.frame(q2), file.path(dir_res, "spls_Q2_by_component.csv"))
# Because the cross-validated Q2 is not positive for the components, no component satisfies the
# Q2 > 0 retention rule; the number of components is therefore FIXED at two for descriptive
# purposes (as stated in the Methods). Q2 per component is reported as obtained.
ncomp <- 2

# 2) keepX / keepY ------------------------------------------------------------------
tune <- mixOmics::tune.spls(X, Y, ncomp = ncomp,
                            test.keepX = seq(5, 50, by = 5), test.keepY = 3:5,
                            validation = "Mfold", folds = NFOLDS, nrepeat = NREPEAT_CV,
                            measure = "cor", mode = "regression",
                            BPPARAM = bp, progressBar = FALSE)
keepX <- tune$choice.keepX
keepY <- tune$choice.keepY
pdf(file.path(dir_fig, "spls_tuning.pdf"), width = 7, height = 5); print(plot(tune)); dev.off()

# 3) final model ---------------------------------------------------------------------
fit <- mixOmics::spls(X, Y, ncomp = ncomp, keepX = keepX, keepY = keepY, mode = "regression")
saveRDS(fit, file.path(dir_res, "spls_model.rds"))
# explained variance of each block per component (Results: X 45%/9%, Y 49%/34%)
write_csv(tibble(component = paste0("Comp", seq_len(ncomp)),
                 X_explained = as.numeric(fit$prop_expl_var$X),
                 Y_explained = as.numeric(fit$prop_expl_var$Y)),
          file.path(dir_res, "spls_explained_variance.csv"))

sel <- bind_rows(lapply(seq_len(ncomp), function(h) {
  v <- mixOmics::selectVar(fit, comp = h)$X$value
  tibble(gene = rownames(v), loading = v$value.var, component = paste0("Comp", h))
}))
write_csv(sel, file.path(dir_res, "spls_selected_genes.csv"))
sel_genes <- unique(sel$gene)

cor_obs <- cor(X[, sel_genes, drop = FALSE], Y)                # gene x trait correlations
write_csv(as.data.frame(cor_obs) |> tibble::rownames_to_column("gene"),
          file.path(dir_res, "spls_gene_trait_correlation.csv"))

# 4) bootstrap: selection frequency and CI of gene-trait correlations ---------------
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED + 2)
boot_one <- function(b) {
  idx <- sample(nrow(X), replace = TRUE)
  m <- try(mixOmics::spls(X[idx, ], Y[idx, ], ncomp = ncomp, keepX = keepX,
                          keepY = keepY, mode = "regression"), silent = TRUE)
  if (inherits(m, "try-error")) return(NULL)
  chosen <- unique(unlist(lapply(seq_len(ncomp), function(h)
    rownames(mixOmics::selectVar(m, comp = h)$X$value))))
  list(chosen = chosen, cor = cor(X[idx, sel_genes, drop = FALSE], Y[idx, ]))
}
boot <- parallel::mclapply(seq_len(BOOT_B), boot_one, mc.cores = N_CORES, mc.set.seed = TRUE)
boot <- Filter(Negate(is.null), boot)
RNGkind("Mersenne-Twister"); set.seed(SEED)

boot_freq <- tibble(gene = sel_genes,
                    boot_frequency = vapply(sel_genes, function(g)



                      mean(vapply(boot, function(r) g %in% r$chosen, logical(1))), numeric(1)))
cor_arr <- simplify2array(lapply(boot, `[[`, "cor"))            # gene x trait x B
ci <- bind_rows(lapply(traits, function(tr) {
  tibble(gene = sel_genes, trait = tr,
         r = cor_obs[, tr],
         ci_low = apply(cor_arr[, tr, ], 1, quantile, 0.025, na.rm = TRUE),
         ci_high = apply(cor_arr[, tr, ], 1, quantile, 0.975, na.rm = TRUE))
})) |> mutate(ci_excludes_zero = ci_low > 0 | ci_high < 0) |>
  left_join(boot_freq, by = "gene")
write_csv(ci, file.path(dir_res, "spls_bootstrap_gene_trait_CI.csv"))

# 5) leave-one-animal-out retention ----------------------------------------------------
loo <- bind_rows(lapply(seq_len(nrow(X)), function(i) {
  m <- mixOmics::spls(X[-i, ], Y[-i, ], ncomp = ncomp, keepX = keepX, keepY = keepY,
                      mode = "regression")
  chosen <- unique(unlist(lapply(seq_len(ncomp), function(h)
    rownames(mixOmics::selectVar(m, comp = h)$X$value))))
  tibble(left_out = rownames(X)[i], gene = sel_genes, retained = sel_genes %in% chosen)
}))
loo_sum <- loo |> group_by(gene) |> summarise(loo_retention = mean(retained), .groups = "drop")
write_csv(loo_sum, file.path(dir_res, "spls_leave_one_out_retention.csv"))

# 6) permutation of the X/Y association (first component) --------------------------------
# Note: keepX/keepY are kept fixed, which is slightly anti-conservative.
assoc <- function(Xm, Ym) {
  m <- mixOmics::spls(Xm, Ym, ncomp = 1, keepX = keepX[1], keepY = keepY[1], mode = "regression")
  abs(cor(m$variates$X[, 1], m$variates$Y[, 1]))
}
obs_assoc <- assoc(X, Y)
RNGkind("L'Ecuyer-CMRG"); set.seed(SEED + 3)
perm_assoc <- unlist(parallel::mclapply(seq_len(N_PERM), function(i) assoc(X, Y[sample(nrow(Y)), ]),
                                        mc.cores = N_CORES, mc.set.seed = TRUE))
RNGkind("Mersenne-Twister"); set.seed(SEED)
p_assoc <- (sum(perm_assoc >= obs_assoc) + 1) / (length(perm_assoc) + 1)
write_csv(tibble(observed_cor_comp1 = obs_assoc, p_permutation = p_assoc, n_perm = N_PERM),
          file.path(dir_res, "spls_permutation_association.csv"))

# 6b) sensitivity analysis: treatment-adjusted carcass traits -----------------------------
# The integrative sPLS pools both treatments, so treatment is a potential confounder. Traits are
# residualized on treatment (z-scored again) and the model refitted with the same keepX/keepY.
trt <- meta$treatment
Y_adj <- scale(apply(Y, 2, function(y) residuals(lm(y ~ trt))))
rownames(Y_adj) <- rownames(X)
fit_adj <- mixOmics::spls(X, Y_adj, ncomp = ncomp, keepX = keepX, keepY = keepY, mode = "regression")
sel_adj <- unique(unlist(lapply(seq_len(ncomp), function(h)
  rownames(mixOmics::selectVar(fit_adj, comp = h)$X$value))))
cor_adj <- cor(X[, sel_genes, drop = FALSE], Y_adj)
write_csv(tibble(n_selected_main = length(sel_genes), n_selected_adjusted = length(sel_adj),
                 n_shared = length(intersect(sel_genes, sel_adj)),
                 prop_main_retained = length(intersect(sel_genes, sel_adj)) / length(sel_genes)),
          file.path(dir_res, "spls_sensitivity_treatment_adjusted_overlap.csv"))
write_csv(as.data.frame(cor_adj) |> tibble::rownames_to_column("gene"),
          file.path(dir_res, "spls_sensitivity_treatment_adjusted_correlation.csv"))

# 7) standard mixOmics plots ---------------------------------------------------------------
pdf(file.path(dir_fig, "spls_correlation_circle.pdf"), width = 6, height = 6)
mixOmics::plotVar(fit, cutoff = 0.70, comp = 1:2, var.names = c(TRUE, TRUE))
dev.off()
pdf(file.path(dir_fig, "spls_indiv.pdf"), width = 7, height = 5)
mixOmics::plotIndiv(fit, ind.names = FALSE, rep.space = "XY-variate")
dev.off()
message("04_spls_carcass.R finished. Report Q2 as obtained (see spls_Q2_by_component.csv).")
