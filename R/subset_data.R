#' Subset The Kinome Data
#'
#' Subset the Kinome data to a given list of peptides and conditions
#'
#' @param chipdata Pamchip-PTK. A Pamchip-PTK Object with all the data available
#' @param peptides character. A vector of peptides of interest
#' @param class character. A specification of the sample class you want to filter to
#'
#' @return A tibble with the filtered and asinh transformed data
#' @export
#'
#' @import dplyr tidyselect
#'
#' @examples
#' TRUE
subset_data <- function(chipdata, peptides, class) {
  if (!class %in% classes(chipdata)) {
    stop("Unrecognized class supplied. Cannot continue.")
  }

  if (!all(peptides %in% peptides(chipdata))) {
    stop("Unrecognized peptides in the peptide list. Cannot continue.")
  }

  samples <- pheno_data(chipdata) %>%
    dplyr::filter(.data$class == {
      {
        class
      }
    }) %>%
    dplyr::pull(.data$sample)

  filtered_data <- exp_data(chipdata) %>%
    dplyr::filter(.data$sample %in% samples) %>%
    dplyr::select(tidyselect::any_of(peptides)) %>%
    dplyr::mutate(dplyr::across(tidyselect::everything(), asinh))

  filtered_data
}
