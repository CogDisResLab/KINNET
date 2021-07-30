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

  post <-
    (probability_dataset[peptide, kinase] * marginal_kinase_probability[kinase]) / marginal_peptide_probability[peptide]

  out <- list(peptide = peptide,
              kinase = kinase,
              posterior = post)

  out
}

#' Generate a Probability Matrix of observing a given peptide/kinase pair
#'
#' @param annotated_df tbl_df. A tibble with a list of all peptide/kinase pairs
#' @param identifier The identifier to use in outputs. Can be either "Gene_Symbol" or "Kinase"
#'
#' @return A tibble with the unadjusted probabilities
#' @import dplyr tidyr purrr
#' @importFrom rlang .data
#'
#' @keywords internal
#' @examples
#' TRUE
generate_probability_matrix <- function(annotated_df, identifier) {
  probability_matrix <- annotated_df %>%
    dplyr::select(.data$ID, .data[[identifier]]) %>%
    dplyr::mutate(value = 1) %>%
    tidyr::pivot_wider(names_from = .data[[identifier]]) %>%
    tibble::column_to_rownames("ID") %>%
    as.matrix()

  probability_matrix[is.na(probability_matrix)] <- 0

  probability <-
    probability_matrix / (reduce(dim(probability_matrix), `*`))

  probability
}

#' Generate Posterior probabilities given priors
#'
#' @param annotated_df tbl_df. A tibble with a list of all peptide/kinase pairs
#' @param identifier The identifier to use in outputs. Can be either "Gene_Symbol" or "Kinase"
#'
#' @return A tibble with posterior probabilities for each peptide/kinase pair
#' @import dplyr purrr
#' @importFrom rlang .data
#'
#' @keywords internal
#' @examples
#' TRUE
generate_posterior_probability_df <-
  function(annotated_df, identifier = "Gene_Symbol") {
    probability_df <-
      generate_probability_matrix(annotated_df, identifier)

    posterior_probs <- annotated_df %>%
      dplyr::select(.data$ID, .data[[identifier]]) %>%
      dplyr::rename(peptide = .data$ID,
                    kinase = .data[[identifier]]) %>%
      purrr::pmap_dfr(~ posterior(..1, ..2, probability_df))

    posterior_probs
  }

#' Calculated Updated Probabilities of each peptide/kinase pair
#'
#' @param chiptype character. Either "STK" or "PTK"
#' @param assignment_df tbl_df. A tibble with the potential kianse assignments
#' @param identifier The identifier to use in outputs. Can be either "Gene_Symbol" or "Kinase"
#'
#' @return an updated dataframe with the probabilities
#' @export
#'
#' @import dplyr
#' @importFrom rlang .data
#' @importFrom utils data
#'
#' @examples
#' TRUE
update_probability_matrix <-
  function(chiptype, assignment_df, identifier = "Gene_Symbol") {
    if (!chiptype %in% c("PTK", "STK")) {
      stop("Incorrect Chip Type specified")
    }

    if (chiptype == "PTK") {
      annotation_posterior <- ptk_probability_matrix
    } else if (chiptype == "STK") {
      annotation_posterior <- stk_probability_matrix
    } else {
      stop("Invalid Chip Specified.\nPlease choose between STK and PTK.")
    }

    if (identifier == "Gene_Symbol") {
      interactome <- kinase_interactome_gene
    } else if (identifier == "Kinase") {
      interactome <- kinase_interactome_kinase
    } else {
      stop("Invalid Identifier Specified.\nPlease choose between Kinase and Gene_Symbol")
    }


    assignment_posterior <-
      generate_posterior_probability_df(assignment_df)

    updated_df <- dplyr::inner_join(
      assignment_posterior,
      annotation_posterior,
      by = c("peptide", "kinase"),
      suffix = c(".assigned", ".reference")
    ) %>%
      mutate(
        foldchange = .data$posterior.assigned / .data$posterior.reference,
        log.assigned = log10(.data$posterior.assigned),
        log.reference = log10(.data$posterior.reference),
        logged.foldchange = .data$log.assigned - .data$log.reference
      )

    updated_df
  }
