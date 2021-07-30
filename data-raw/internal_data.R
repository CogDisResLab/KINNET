## code to prepare `internal_data` dataset goes here

source("data-raw/ptk_annotation.R")
source("data-raw/stk_annotation.R")
source("data-raw/kinase_interactome_gene.R")
source("data-raw/kinase_interactome_kinase.R")
source("data-raw/ptk_probability_matrix.R")
source("data-raw/stk_probability_matrix.R")

usethis::use_data(kinase_interactome_gene,
                  kinase_interactome_kinase,
                  ptk_annotation,
                  stk_annotation,
                  ptk_probability_matrix_gene,
                  ptk_probability_matrix_kinase,
                  stk_probability_matrix_gene,
                  stk_probability_matrix_kinase,
                  internal = TRUE, overwrite = TRUE)
