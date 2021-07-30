## code to prepare `ptk_probability_matrix` dataset goes here

ptk_probability_matrix_gene <- KINNET::generate_posterior_probability_df(ptk_annotation)
ptk_probability_matrix_kinase <- KINNET::generate_posterior_probability_df(ptk_annotation, "Kinase")
