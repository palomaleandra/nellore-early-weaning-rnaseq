# ==============================================================================
# 06_enrichment_overlap.R
# Over-representation analysis (enrichR) of sPLS-DA / sPLS gene sets and overlap
# between models (Euler diagram, Fig. 5). ORA on pre-selected genes is contextual
# evidence only.
# Inputs : results/splsda_<tp>_selected.csv, results/spls_selected_genes.csv
# Outputs: results/ora_*.csv, results/overlap_*.csv, figures/euler_overlap.pdf
# Requires internet access (Enrichr web service).
# ==============================================================================
source(here::here("R", "00_setup.R"))
suppressPackageStartupMessages({ library(enrichR); library(eulerr); library(dplyr); library(readr) })

DBS <- resolve_enrichr_dbs()          # newest GO-BP / KEGG / Reactome libraries available
writeLines(DBS, file.path(dir_res, "enrichr_libraries_used.txt"))   # report in the Methods

genes_of <- function(file, comp = NULL) {
  d <- read_csv(file.path(dir_res, file), show_col_types = FALSE)
  if (!is.null(comp)) d <- d |> filter(component %in% comp)
  unique(d$gene)
}
wean <- genes_of("splsda_weaning_selected.csv")
slaug <- genes_of("splsda_preslaughter_selected.csv")
carc <- genes_of("spls_selected_genes.csv")

ora <- function(genes, label) {
  if (length(genes) < 2) return(invisible(NULL))
  res <- enrichR::enrichr(genes, DBS)
  out <- bind_rows(lapply(names(res), function(db) {
    r <- res[[db]]
    if (is.null(r) || nrow(r) == 0) return(NULL)
    r$database <- db
    r
  })) |>
    filter(P.value < 0.05) |>                              # nominal P < 0.05 (exploratory), as in the Methods
    mutate(significant_FDR = Adjusted.P.value < 0.05) |>   # BH q-value within each database
    arrange(P.value)
  write_csv(out, file.path(dir_res, sprintf("ora_%s.csv", label)))
}

ora(wean, "weaning_splsda")
ora(slaug, "preslaughter_splsda")
ora(carc, "spls_carcass")

# ---- overlaps --------------------------------------------------------------------------------
ov <- list(
  `sPLS-DA weaning`      = wean,
  `sPLS-DA pre-slaughter` = slaug,
  `sPLS carcass`         = carc
)
inter <- list(
  weaning_x_preslaughter = intersect(wean, slaug),
  preslaughter_x_carcass = intersect(slaug, carc),
  weaning_x_carcass       = intersect(wean, carc),
  all_three               = Reduce(intersect, ov)
)
for (n in names(inter)) {
  write_csv(tibble::tibble(gene = inter[[n]]), file.path(dir_res, sprintf("overlap_%s.csv", n)))
}
ora(inter$weaning_x_preslaughter, "overlap_weaning_x_preslaughter")
ora(inter$preslaughter_x_carcass, "overlap_preslaughter_x_carcass")

pdf(file.path(dir_fig, "euler_overlap.pdf"), width = 5, height = 4)
print(plot(eulerr::euler(ov), quantities = TRUE, fills = list(alpha = 0.6)))
dev.off()
message("06_enrichment_overlap.R finished.")
