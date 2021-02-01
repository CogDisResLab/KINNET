## code to prepare `selected_peptides` dataset goes here

library(tidyverse)

selected_list <- readRDS("~/experiments/justin-kinome-modelling/kinome_data/changed_peps_forEachComparison_JustinPancKinome.RDS")

# PANC1 Comparisons
panc1_wt <- tibble(cell_line = "PANC1", reference = "PANC1", comparison = "WT", peptides = selected_list$panc1vsBenign)
panc1_dig <- tibble(cell_line = "PANC1", reference = "PANC1", comparison = "Dig", peptides = selected_list$panc1vsDig)
panc1_met <- tibble(cell_line = "PANC1", reference = "PANC1", comparison = "Met", peptides = selected_list$panc1vsMet)
panc1_sim <- tibble(cell_line = "PANC1", reference = "PANC1", comparison = "Sim", peptides = selected_list$panc1vsSim)

# PDCL5 Comparisons
pdcl5_wt <- tibble(cell_line = "PDCL5", reference = "PDCL5", comparison = "WT", peptides = selected_list$pdcl5vsBenign)
pdcl5_dig <- tibble(cell_line = "PDCL5", reference = "PDCL5", comparison = "Dig", peptides = selected_list$pdcl5vsDig)
pdcl5_met <- tibble(cell_line = "PDCL5", reference = "PDCL5", comparison = "Met", peptides = selected_list$pdcl5vsMet)
pdcl5_sim <- tibble(cell_line = "PDCL5", reference = "PDCL5", comparison = "Sim", peptides = selected_list$pdcl5vsSim)

# PDCL15 Comparisons
pdcl15_wt <- tibble(cell_line = "PDCL15", reference = "PDCL15", comparison = "WT", peptides = selected_list$pdcl15vsBenign)
pdcl15_dig <- tibble(cell_line = "PDCL15", reference = "PDCL15", comparison = "Dig", peptides = selected_list$pdcl15vsDig)
pdcl15_met <- tibble(cell_line = "PDCL15", reference = "PDCL15", comparison = "Met", peptides = selected_list$pdcl15vsMet)
pdcl15_sim <- tibble(cell_line = "PDCL15", reference = "PDCL15", comparison = "Sim", peptides = selected_list$pdcl15vsSim)

selected_peptides <- bind_rows(
  panc1_wt,
  panc1_dig,
  panc1_met,
  panc1_sim,
  pdcl5_wt,
  pdcl5_dig,
  pdcl5_met,
  pdcl5_sim,
  pdcl15_wt,
  pdcl15_dig,
  pdcl15_met,
  pdcl15_sim
)

usethis::use_data(selected_peptides, overwrite = TRUE)
