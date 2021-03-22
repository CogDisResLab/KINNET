#' Calculate the Structural Hamming Distance between two kinased graphs
#'
#' @param reference a reference result from assign_kinases()
#' @param comparison a comparison result from assign_kinases()
#'
#' @return A positive integer
#' @export
#'
#' @importFrom bnlearn shd
#'
#' @examples
#' TRUE
kinnet_shd <- function(reference, comparison) {

  pair <- equalize_kinase_graphs(reference, comparison)

  res <- bnlearn::shd(pair$reference, pair$comparison)

  res
}

#' Calculate the Hamming Distance between two kinased graphs
#'
#' @param reference a reference result from assign_kinases()
#' @param comparison a comparison result from assign_kinases()
#'
#' @return A positive integer
#' @export
#'
#' @importFrom bnlearn hamming
#'
#' @examples
#' TRUE
kinnet_hamming <- function(reference, comparison) {

  pair <- equalize_kinase_graphs(reference, comparison)

  res <- bnlearn::hamming(pair$reference, pair$comparison)

  res
}


#' Calculate Normalized Structural Hamming Distance
#'
#' @param reference a reference result from assign_kinases()
#' @param comparison a comparison result from assign_kinases()
#'
#' @return A positive integer
#' @export
#'
#' @importFrom bnlearn shd
#'
#' @examples
#' TRUE
kinnet_normalized_shd <- function(reference, comparison) {

  pair <- equalize_kinase_graphs(reference, comparison)

  shd <- bnlearn::shd(pair$reference, pair$comparison)

  arcnum <- map_dbl(pair, narcs) %>% sum

  norm_shd <-round(shd/arcnum, 5)

  norm_shd
}


#' Calculate Normalized Hamming Distance
#'
#' @param reference a reference result from assign_kinases()
#' @param comparison a reference result from assign_kinases()
#'
#' @return A positive integer
#' @export
#'
#' @importFrom bnlearn hamming
#'
#' @examples
#' TRUE
kinnet_normalized_hamming <- function(reference, comparison) {

  pair <- equalize_kinase_graphs(reference, comparison)

  hamming <- bnlearn::shd(pair$reference, pair$comparison)

  arcnum <- map_dbl(pair, narcs) %>% sum

  norm_hamming <-round(shd/arcnum, 5)

  norm_hamming
}

