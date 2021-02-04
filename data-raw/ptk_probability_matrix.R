## code to prepare `ptk_probability_matrix` dataset goes here

data("ptk_annotation")

ptk_probability_matrix <- KINNET::generate_posterior_probability_df(ptk_annotation)
