## code to prepare `internal_data` dataset goes here

source("data-raw/kinase_interactome.R")
source("data-raw/ptk_annotation.R")
source("data-raw/stk_annotation.R")
source("data-raw/ptk_probability_matrix.R")
source("data-raw/stk_probability_matrix.R")

usethis::use_data(kinase_interactome,
                  ptk_annotation,
                  stk_annotation,
                  ptk_probability_matrix,
                  stk_probability_matrix,
                  internal = TRUE, overwrite = TRUE)
