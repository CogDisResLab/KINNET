## code to prepare `kinase_interactome` dataset goes here
##
library(tidyverse)
library(devtools)

interactome_file <- file.path("data-raw/data/RegPhos_kinase_PPI_human.txt")

interactome <- read_tsv(interactome_file) %>%
  select(GENE_a, GENE_b) %>%
  rename(from = GENE_a,
         to = GENE_b)

interactome_na_to <- interactome$to %>%
  tibble(from = NA, to = .)

interactome_na_from <- interactome$from %>%
  tibble(from = ., to = NA)

kinase_interactome <- bind_rows(interactome, interactome_na_to, interactome_na_from)

usethis::use_data(kinase_interactome, overwrite = TRUE)
