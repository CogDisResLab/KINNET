## code to prepare `stk_probability_matrix` dataset goes here

data("stk_annotation")

stk_probability_matrix <- KINNET::generate_posterior_probability_df(stk_annotation)
