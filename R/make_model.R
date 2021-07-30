#' Generate a bayes net model
#'
#' @param expression tbl_df. A tibble with the activity data from the kinome chip.
#' @param iterations numeric. Number of iterations to run the mode.
#' @param threshold numeric. Threshold to use for averaging the network
#' @param cluster cluster. (Optional) a cluster from the package parallel
#'
#' @return A list with the strength network dataframe, an averaged network and the threshold used to generate that averaged network.
#' @export
#'
#' @import bnlearn
#'
#' @examples
#' TRUE
make_model <-
  function(expression,
           iterations = 200,
           threshold = NULL,
           cluster = NULL) {
    strength_model <- bnlearn::boot.strength(expression,
                                             R = iterations,
                                             algorithm = "hc",
                                             cluster = cluster)

    significance_threshold <- if_else(is.null(threshold),
                                      attr(strength_model, "threshold"),
                                      threshold)

    out <- list(
      strength_net = strength_model,
      avg_net = bnlearn::averaged.network(strength_model, threshold = significance_threshold),
      threshold = significance_threshold
    )

    out
  }
