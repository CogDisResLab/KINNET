#' PTK Peptide to Kinase Probability Matrix
#'
#' A dataaframe with the bayesian posteriors of finding a given peptide/kinase pair
#'
#' @format A dataframe with 991 rows and 3 columns
#'
#' \describe{
#' \item{peptide}{The peptide ID}
#' \item{kinase}{Gene Symbol of the kinase}
#' \item{posterior}{Posterior probability of the pair}
#' }
#'
"ptk_probability_matrix"
