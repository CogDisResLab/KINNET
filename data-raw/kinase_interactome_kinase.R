## code to prepare `kinase_interactome_kinase` dataset goes here
##
library(tidyverse)
library(devtools)

interactome_file <-
  file.path("data-raw", "data", "RegPhos_kinase_PPI_human.txt")
stk_annotation_file <-
  file.path("data-raw", "data", "STK_Annotation.csv")
ptk_annotation_file <-
  file.path("data-raw", "data", "PTK_Annotation.csv")

interactome <- read_tsv(interactome_file, col_types = cols(.default = col_character())) |>
  select(GENE_a, GENE_b) |>
  rename(from = GENE_a,
         to = GENE_b)

stk_annotation_loaded <-
  read_csv(stk_annotation_file, col_types = cols(.default = col_character())) |>
  select(Kinase, Gene_Symbol)

stk_interactome_kinase <- interactome |>
  inner_join(stk_annotation_loaded, by = c(from = "Gene_Symbol")) |>
  select(-from) |>
  rename(from = Kinase) |>
  inner_join(stk_annotation_loaded, by = c(to = "Gene_Symbol")) |>
  select(-to) |>
  rename(to = Kinase) |>
  unique()

ptk_annotation_loaded <-
  read_csv(ptk_annotation_file, col_types = cols(.default = col_character())) |>
  select(Kinase, Gene_Symbol)

ptk_interactome_kinase <- interactome |>
  inner_join(ptk_annotation_loaded, by = c(from = "Gene_Symbol")) |>
  select(-from) |>
  rename(from = Kinase) |>
  inner_join(ptk_annotation_loaded, by = c(to = "Gene_Symbol")) |>
  select(-to) |>
  rename(to = Kinase) |>
  unique()

full_interactome <-
  bind_rows(stk_interactome_kinase, ptk_interactome_kinase)

interactome_na_to <- full_interactome$to |>
  (\(x) tibble(from = NA, to = x))()

interactome_na_from <- full_interactome$from |>
  (\(x) tibble(from = x, to = NA))()

kinase_interactome_kinase <-
  bind_rows(full_interactome, interactome_na_to, interactome_na_from)
