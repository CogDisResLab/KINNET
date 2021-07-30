## code to prepare `stk_probability_matrix` dataset goes here

stk_probability_matrix_gene <- KINNET::generate_posterior_probability_df(stk_annotation)
stk_probability_matrix_kinase <- KINNET::generate_posterior_probability_df(stk_annotation, "Kinase")
