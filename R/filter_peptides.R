#' Do a standardized fit of expression data from a single peptide's activity
#'
#' @param expr a dataframe containing the activity data of only one peptide
#'
#' @return a list of three items std_10, std_50 and std_200, of standardized linear coefficients
#' @export
#'
#' @importFrom dplyr filter group_by group_split
#' @importFrom effectsize standardize_parameters
#' @importFrom purrr map map_dbl
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
fit_standardized <- function(expr) {
  res <- expr %>%
    dplyr::filter(.data$exposure %in% c(10, 50, 200)) %>%
    dplyr::group_by(.data$exposure) %>%
    dplyr::group_split() %>%
    purrr::map(~ lm(activity ~ cycle, data = .x)) %>%
    purrr::map(~ effectsize::standardize_parameters(.x)) %>%
    purrr::map("Std_Coefficient") %>%
    purrr::map(2) %>%
    purrr::map_dbl(~ round(.x, 3))

  names(res) <- c("std_10", "std_50", "std_200")
  output <- as.list(res)
  output
}

#' Return a list of significant peptides
#'
#' @param chipdata an object of class PamchipSTK or PamchipPTK
#' @param threshold a lower cutoff for significance
#'
#' @return a list of significant peptides for each class
#' @export
#'
#' @importFrom dplyr mutate across inner_join ungroup group_by filter summarise group_split if_else
#' @importFrom purrr map
#' @importFrom tidyr nest unnest_wider pivot_longer
#' @importFrom tidyselect vars_select_helpers
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
filter_peptides <- function(chipdata, threshold = 0.75) {
  pheno <- pheno_data(chipdata) %>%
    dplyr::mutate(exposure = as.numeric(.data$exposure),
           cycle = as.numeric(.data$cycle) - 92)

  expr <- exp_data(chipdata) %>%
    dplyr::mutate(dplyr::across(tidyselect::vars_select_helpers$where(is.numeric), ~ dplyr::if_else(.x < 0, 0, .x)))

  nclasses <- length(classes(chipdata))
  nreplicates <- round(12/nclasses)


  augmented <-
    pheno %>%
    dplyr::inner_join(expr, by = "sample") %>%
    tidyr::pivot_longer(matches("\\_\\s?\\d+\\_", perl = TRUE), names_to = "peptide", values_to = "activity") %>%
    dplyr::group_by(.data$barcode, .data$array, .data$class, .data$peptide) %>%
    nest() %>%
    dplyr::mutate(cors = purrr::map(.data$data, fit_standardized)) %>%
    tidyr::unnest_wider(col = c(.data$cors)) %>%
    dplyr::ungroup()

  filtered <- augmented %>%
    dplyr::group_by(.data$class, .data$peptide) %>%
    dplyr::filter(.data$std_10 > threshold, .data$std_50 > threshold, .data$std_200 > threshold) %>%
    dplyr::summarise(n = n()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(.data$n == nreplicates) %>%
    dplyr::group_by(class) %>%
    dplyr::group_split() %>%
    purrr::map(~ dplyr::pull(.x, .data$peptide)) %>%
    purrr::map(~ unique(.x))

  names(filtered) <- classes(chipdata)

  filtered
}
