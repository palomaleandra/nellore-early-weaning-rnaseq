# ==============================================================================
# 01_phenotypes.R
# Performance and carcass traits (Table 1).
#    - body weight: linear mixed model, repeated measures
#    - HCW, CY, REA: linear mixed models (REA fitted on log scale)
#    - SFT, intramuscular lipids: Gamma GLMM, log link
#    - EMMs are reported on the ORIGINAL scale (type = "response")
#    - Benjamini-Hochberg adjustment across the five carcass traits
#    - Monte Carlo power (simr) and TOST equivalence of the RNA-seq subset
# Inputs : data/metadata/phenotypes.csv, data/metadata/bodyweight_long.csv
# Outputs: results/table1_carcass.csv, results/table1_adg.csv, results/table_S5_tost.csv, table_S6_effects_power.csv
# ==============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({
   library(lme4); library(lmerTest); library(emmeans); library(dplyr); library(readr)
})

# phenotypes.csv columns (one row per animal):
#   animal_id, treatment (CW/EW), paddock, rnaseq (TRUE/FALSE),
#   adg_120_205, adg_205_914, adg_120_914, bw_120, bw_205, bw_914, hcw, cy, rea, sft, lipids
pheno <- load_phenotypes()

# ---- helper: EMMs back-transformed + contrast + Cohen's d ---------------------
summarise_model <- function(model, trait, data) {
  emm <- emmeans(model, ~ treatment, type = "response")
  est <- as.data.frame(emm)
  est_col <- intersect(c("emmean", "response"), names(est))[1]
  ctr <- as.data.frame(pairs(emm, reverse = TRUE, type = "response")) # EW - CW
  eff <- as.data.frame(eff_size(emm, sigma = sigma(model), edf = nrow(data) - 2))
  # eff_size() contrasts CW - EW (factor-level order): same sign convention as Supplementary
  # Table S6 (negative d = CW lower than EW).
  tibble(
    trait    = trait,
    CW_emm   = est[[est_col]][est$treatment == "CW"],
    CW_se    = est$SE[est$treatment == "CW"],
    EW_emm   = est[[est_col]][est$treatment == "EW"],
    EW_se    = est$SE[est$treatment == "EW"],
    contrast = ctr[[2]],            # difference (identity link) or ratio (log link)
    contrast_se = ctr$SE,
    df       = ctr$df,
    stat     = ctr[[intersect(c("t.ratio", "z.ratio"), names(ctr))[1]]],   # t (LMM) or z (GLMM)
    p_value = ctr$p.value,
    cohen_d = eff$effect.size[1],
    cohen_d_low = eff$lower.CL[1],
    cohen_d_high = eff$upper.CL[1],
    AIC      = AIC(model)
  )
}

# ---- carcass traits -----------------------------------------------------------
m_hcw <- lmer(hcw ~ treatment + (1 | paddock), data = pheno)
m_cy   <- lmer(cy ~ treatment + (1 | paddock), data = pheno)
m_rea <- lmer(log(rea) ~ treatment + (1 | paddock), data = pheno)    # back-transformed by emmeans
m_sft <- glmer(sft     ~ treatment + (1 | paddock), data = pheno, family = Gamma(link = "log"))
m_lip <- glmer(lipids ~ treatment + (1 | paddock), data = pheno, family = Gamma(link = "log"))
# Note: with few paddocks the random effect can be singular (isSingular(m_*)).

carcass <- bind_rows(
  summarise_model(m_hcw, "HCW (kg)",   pheno),
  summarise_model(m_cy, "CY (%)",      pheno),
  summarise_model(m_rea, "REA (cm2)", pheno),
  summarise_model(m_sft, "SFT (mm)",   pheno),
  summarise_model(m_lip, "Lipids (%)", pheno)
) |>
  mutate(p_BH = p.adjust(p_value, method = "BH"))   # multiplicity across carcass traits
write_csv(carcass, file.path(dir_res, "table1_carcass.csv"))

# ---- average daily gain -------------------------------------------------------
adg_traits <- c("adg_205_914", "adg_120_205", "adg_120_914")
adg_models <- setNames(lapply(adg_traits, function(v)
  lm(reformulate("treatment", response = v), data = pheno)), adg_traits)
adg <- bind_rows(Map(function(m, v) summarise_model(m, v, pheno), adg_models, adg_traits))
write_csv(adg, file.path(dir_res, "table1_adg.csv"))

# ---- body weight, repeated measures ------------------------------------------
bw <- readr::read_csv(file.path(dir_meta, "bodyweight_long.csv"), show_col_types = FALSE)
bw$treatment <- factor(bw$treatment, levels = c("CW", "EW"))
bw$age <- factor(bw$age_days)
# Homogeneous vs heterogeneous residual variance (varIdent by age), compared by AIC (REML fits
# with identical fixed effects). nlme is used because lme4 cannot model residual heteroscedasticity.
m_bw_hom <- nlme::lme(bw ~ treatment * age, random = ~ 1 | animal_id, data = bw, method = "REML")
m_bw_het <- nlme::lme(bw ~ treatment * age, random = ~ 1 | animal_id, data = bw, method = "REML",
                      weights = nlme::varIdent(form = ~ 1 | age))
aic_bw <- AIC(m_bw_hom, m_bw_het)
shap <- shapiro.test(residuals(m_bw_het, type = "normalized"))
write.csv(data.frame(model = c("Homogeneous variance", "Heterogeneous variance"), AIC = aic_bw$AIC,
                     Shapiro_W = shap$statistic, Shapiro_p = shap$p.value),



          file.path(dir_res, "bodyweight_AIC_variance_structures.csv"), row.names = FALSE)   # Table S4
m_bw <- if (AIC(m_bw_het) < AIC(m_bw_hom)) m_bw_het else m_bw_hom
capture.output(anova(m_bw), file = file.path(dir_res, "bodyweight_anova.txt"))
write_csv(as.data.frame(pairs(emmeans(m_bw, ~ treatment | age), reverse = TRUE)),
          file.path(dir_res, "bodyweight_contrasts.csv"))

# ---- marginal and conditional R2 (Nakagawa & Schielzeth 2013) ------------------
r2_tab <- bind_rows(lapply(list(HCW = m_hcw, CY = m_cy, REA = m_rea, SFT = m_sft, Lipids = m_lip),
                           function(m) { r <- performance::r2_nakagawa(m)
                             tibble(R2_marginal = as.numeric(r$R2_marginal),
                                    R2_conditional = as.numeric(r$R2_conditional)) }),
                    .id = "trait")
write_csv(r2_tab, file.path(dir_res, "table_R2_nakagawa.csv"))

# ---- Monte Carlo power (simr), all five carcass traits --------------------------
# Power to detect the observed treatment effect, 1,000 simulations.
power_one <- function(model) {
  tst <- if (inherits(model, "glmerMod")) "z" else "t"
  simr::powerSim(model, test = simr::fixed("treatmentEW", tst), nsim = 1000, progress = FALSE)
}
set.seed(SEED)
pw <- lapply(c(list(HCW = m_hcw, CY = m_cy, REA = m_rea, SFT = m_sft, Lipids = m_lip), adg_models),
             power_one)
capture.output(for (n in names(pw)) { cat("
==", n, "==
"); print(pw[[n]]) },
               file = file.path(dir_res, "power_simr.txt"))

# ---- Supplementary Table S6: AIC, Cohen's d (95% CI) and power ----------------------
power_tab <- tibble(trait = names(pw), power = vapply(pw, function(p) summary(p)$mean, numeric(1)))
s6 <- bind_rows(carcass, adg) |> select(trait, AIC, cohen_d, cohen_d_low, cohen_d_high) |>
  mutate(trait = recode(trait, "HCW (kg)" = "HCW", "CY (%)" = "CY", "REA (cm2)" = "REA",
                        "SFT (mm)" = "BFT", "Lipids (%)" = "IMF")) |>
  left_join(power_tab |> mutate(trait = recode(trait, SFT = "BFT", Lipids = "IMF")), by = "trait")
write_csv(s6, file.path(dir_res, "table_S6_effects_power.csv"))

# ---- TOST: is the RNA-seq subsample equivalent to the remaining animals? -------
# Supplementary Table S5: selected (n = 16, rnaseq = TRUE) vs remaining animals (n = 24), all 40
# animals pooled across treatments. Equivalence bounds = +/- 1 SD of the trait in all animals;
# p_equiv = larger of the two one-sided p-values (Welch). Includes body weights, ADG, carcass
# traits and IMF.
tost_traits <- c(bw_120 = "BW at 120 d", bw_205 = "BW at 205 d", bw_914 = "BW at 914 d",
                 adg_205_914 = "ADG (205-914)", adg_120_205 = "ADG (120-205)",
                 adg_120_914 = "ADG (120-914)", hcw = "HCW", cy = "CY", rea = "REA",
                 sft = "BFT", lipids = "IMF")
tost <- bind_rows(lapply(names(tost_traits), function(v) {
  sel <- pheno[[v]][pheno$rnaseq]; rem <- pheno[[v]][!pheno$rnaseq]
  delta <- sd(pheno[[v]], na.rm = TRUE)
  r <- TOSTER::t_TOST(x = sel, y = rem, eqb = delta, hypothesis = "EQU", var.equal = FALSE)
  tibble(outcome = tost_traits[[v]], selected_mean = mean(sel, na.rm = TRUE),
         remaining_mean = mean(rem, na.rm = TRUE), delta = delta,
         p_equiv = max(r$TOST$p.value[-1]))
}))
# Multiplicity: the claim is that the subsample is equivalent for ALL outcomes, a conjunction
# tested by the intersection-union principle: each TOST is run at alpha = 0.05 without adjustment
# and the overall P is the maximum of the individual P values.
tost$p_equiv_overall_IUT <- max(tost$p_equiv)
write_csv(tost, file.path(dir_res, "table_S5_tost.csv"))
message("01_phenotypes.R finished.")
