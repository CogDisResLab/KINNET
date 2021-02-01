## code to prepare `stk_probability_matrix` dataset goes here

data("stk_annotation")

stk_probability_matrix <- JustinKinomeModelling::generate_posterior_probability_df(stk_annotation)

usethis::use_data(stk_probability_matrix, overwrite = TRUE)
