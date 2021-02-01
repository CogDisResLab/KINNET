#' List of significantly differentially expressed peptides in the dataset
#'
#' A list of all the peptides found to be differentially phosphorylated
#' in each pair of samples
#'
#' @format A dataframe with 833 rows and 4 columns
#'
#' \describe{
#' \item{cell_line}{The cell line for the experiment}
#' \item{reference}{The reference sample}
#' \item{comparison}{The comparison sample}
#' \item{peptides}{The peptides selected for that pairing}
#' }
#'
"selected_peptides"
