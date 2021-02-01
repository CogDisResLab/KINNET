#' Generate Posterior Probability of seeing a given kinase-peptide pair
#'
#' @param peptide character. ID of the peptide
#' @param kinase character. Gene Symbol of the Kinase
#' @param probability_dataset tbl_df. The aadjacency matrix dataframe with all kinase-peptide pairs
#'
#' @return a list with the peptide, the kinase and the posterior probability
#' @keywords internal
#'
#' @examples
#' TRUE
posterior <- function(peptide, kinase, probability_dataset) {

  marginal_kinase_probability <- colSums(probability_dataset)
  marginal_peptide_probability <- rowSums(probability_dataset)

  post <- (probability_dataset[peptide, kinase] * marginal_kinase_probability[kinase]) / marginal_peptide_probability[peptide]

  out <- list(peptide = peptide,
              kinase = kinase,
              posterior = post)

  out
}

#' Generate a Probability Matrix of observing a given peptide/kinase pair
#'
#' @param annotated_df tbl_df. A tibble with a list of all peptide/kinase pairs
#'
#' @return A tibble with the unadjusted probabilities
#' @import dplyr tidyr purrr
#' @importFrom rlang .data
#'
#' @keywords internal
#' @examples
#' TRUE
generate_probability_matrix <- function(annotated_df) {
  probability_matrix <- annotated_df %>%
    dplyr::select(.data$ID, .data$Gene_Symbol) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = .data$Gene_Symbol) %>%
    tibble::column_to_rownames("ID") %>%
    as.matrix()

  probability_matrix[is.na(probability_matrix)] <- 0

  probability <- probability_matrix/(reduce(dim(probability_matrix), `*`))

  probability
}

#' Generate Posterior probabilities given priors
#'
#' @param annotated_df tbl_df. A tibble with a list of all peptide/kinase pairs
#'
#' @return A tibble with posterior probabilities for each peptide/kinase pair
#' @import dplyr purrr
#' @importFrom rlang .data
#'
#' @keywords internal
#' @examples
#' TRUE
generate_posterior_probability_df <- function(annotated_df) {

  probability_df <- generate_probability_matrix(annotated_df)

  posterior_probs <- annotated_df %>%
    dplyr::select(.data$ID, .data$Gene_Symbol) %>%
    dplyr::rename(peptide = .data$ID,
                  kinase = .data$Gene_Symbol) %>%
    purrr::pmap_dfr(~ posterior(..1, ..2, probability_df))

  posterior_probs
}

#' Calculated Updated Probabilities of each peptide/kinase pair
#'
#' @param chiptype character. Either "STK" or "PTK"
#' @param assignment_df tbl_df. A tibble with the potential kianse assignments
#'
#' @return an updated dataframe with the probabilities
#' @export
#'
#' @import dplyr
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
update_probability_matrix <- function(chiptype, assignment_df) {

  if (!chiptype %in% c("PTK", "STK")) {
    stop("Incorrect Chip Type specified")
  }

  if (chiptype == "PTK") {
    annotation_posterior <- JustinKinomeModelling::ptk_probability_matrix
  } else if (chiptype == "STK") {
    annotation_posterior <- JustinKinomeModelling::stk_probability_matrix
  }


  assignment_posterior <- generate_posterior_probability_df(assignment_df)

  updated_df <- dplyr::inner_join(assignment_posterior,
                                  annotation_posterior,
                                  by = c("peptide", "kinase"),
                                  suffix = c(".assigned", ".reference")) %>%
    mutate(
      foldchange = .data$posterior.assigned/.data$posterior.reference,
      log.assigned = log10(.data$posterior.assigned),
      log.reference = log10(.data$posterior.reference),
      logged.foldchange = .data$log.assigned - .data$log.reference
    )

  updated_df
}
