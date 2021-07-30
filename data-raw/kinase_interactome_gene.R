## code to prepare `kinase_interactome_gene` dataset goes here
##
library(tidyverse)
library(devtools)

interactome_file <- file.path("data-raw", "data", "RegPhos_kinase_PPI_human.txt")

interactome <- read_tsv(interactome_file, col_types = cols(.default = col_character())) |>
  select(GENE_a, GENE_b)  |>
  rename(from = GENE_a,
         to = GENE_b)

interactome_na_to <- interactome$to |>
  (\(x) tibble(from = NA, to = x))()

interactome_na_from <- interactome$from |>
  (\(x) tibble(from = x, to = NA))()

kinase_interactome_gene <- bind_rows(interactome, interactome_na_to, interactome_na_from)
