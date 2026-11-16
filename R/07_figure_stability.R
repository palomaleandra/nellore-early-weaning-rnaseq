# ==============================================================================
# 07_figure_stability.R
# Selection-frequency (stability) figure for the sPLS-DA models (selection-frequency plots).
# Bars = genes of the FINAL model, coloured by direction (mean VST in EW minus CW),
# dashed line = stability threshold.
# Inputs : results/splsda_<tp>_selected.csv
# Outputs: figures/stability_<tp>.pdf / .png
# ==============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(dplyr); library(readr); library(ggplot2) })

plot_stability <- function(tp) {
  d <- read_csv(file.path(dir_res, sprintf("splsda_%s_selected.csv", tp)), show_col_types = FALSE) |>
    mutate(frequency = tidyr::replace_na(frequency, 0)) |>   # 'direction' comes from group means (script 03)
    group_by(component) |>
    arrange(desc(frequency), .by_group = TRUE) |>
    mutate(rank = row_number()) |>
    ungroup()

    lab <- d |> filter(frequency >= STAB_THRESHOLD)
    ggplot(d, aes(rank, frequency, fill = direction)) +
      geom_col(width = 0.85) +
      geom_hline(yintercept = STAB_THRESHOLD, linetype = "dashed", colour = "firebrick") +
      geom_text(data = lab, aes(label = gene), vjust = -0.4, size = 2.6, fontface = "bold") +
      facet_wrap(~ component, scales = "free", nrow = 1) +
      scale_fill_manual(values = c("Higher in EW" = "#E07B54", "Higher in CW" = "#5B9BD5"), name = NULL) +
      labs(x = "Genes of the final model, ranked by selection frequency",
           y = "Selection frequency (f)") +
      theme_classic(base_size = 10) +
      theme(legend.position = "bottom", axis.text.x = element_blank(), axis.ticks.x = element_blank())
}

for (tp in c("weaning", "preslaughter")) {
  p <- plot_stability(tp)
  ggsave(file.path(dir_fig, sprintf("stability_%s.pdf", tp)), p, width = 8, height = 3.6)
  ggsave(file.path(dir_fig, sprintf("stability_%s.png", tp)), p, width = 8, height = 3.6, dpi = 300)
}
message("07_figure_stability.R finished.")
