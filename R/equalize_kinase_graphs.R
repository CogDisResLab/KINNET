#' Equalize the nodesets of two kinased graphs
#'
#' @param reference a reference result from assign_kinases()
#' @param comparison a comparison result from assign_kinases()
#'
#' @return A list with two elements, reference and comparison, that have the
#' same nodesets
#' @export
#'
#' @importFrom bnlearn nodes
#'
#' @examples
#' TRUE
equalize_kinase_graphs <- function(reference, comparison) {
  ref_kinased <-
    render_reduced_kinased_graph(reference, title = "Reference")
  ref_bnnet <- ref_kinased$bnnet
  comp_kinased <-
    render_reduced_kinased_graph(comparison, title = "Comparison")
  comp_bnnet <- comp_kinased$bnnet

  all_nodes <-
    unique(c(bnlearn::nodes(ref_bnnet), bnlearn::nodes(comp_bnnet)))

  ref_bnnet_augmented <- ref_bnnet
  for (node in all_nodes) {
    ref_bnnet_augmented <- safely_add_node(ref_bnnet_augmented, node)
  }

  comp_bnnet_augmented <- comp_bnnet
  for (node in all_nodes) {
    comp_bnnet_augmented <- safely_add_node(comp_bnnet_augmented, node)
  }

  out <- list(reference = ref_bnnet_augmented,
              comparison = comp_bnnet_augmented)

  out
}
