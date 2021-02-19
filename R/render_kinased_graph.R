#' Generate a kianse graph
#'
#' @param analysis_result Object output from assign_kinases()
#' @param title Title of the graph
#' @param render logical. Whether a graph should be rendered or not.
#'
#' @return a graph object
#' @export
#'
#' @import bnlearn dplyr
#' @importFrom rlang .data
#' @importFrom graph graphRenderInfo nodeRenderInfo edgeRenderInfo
#' @importFrom Rgraphviz layoutGraph renderGraph
#'
#' @examples
#' TRUE
render_kinased_graph <- function(analysis_result, title, render = FALSE) {
  mapped_kinases <- analysis_result$mapped_kinases %>%
    dplyr::select(.data$peptide, .data$kinase) %>%
    dplyr::group_by(.data$peptide) %>%
    dplyr::slice_head(n = 1) %>%
    tibble::deframe()

  network_graph <- bnlearn::graphviz.plot(analysis_result$net,
                                          shape = "ellipse", layout = "fdp",
                                          main = title)


  gr <- Rgraphviz::layoutGraph(network_graph, attrs = list(graph = list(rankdir = "TB")))
  graph::graphRenderInfo(gr)$fontsize <- 24
  graph::nodeRenderInfo(gr)$label <- mapped_kinases
  graph::nodeRenderInfo(gr)$fontsize <- 14
  graph::nodeRenderInfo(gr)$fixedsize <- TRUE
  graph::nodeRenderInfo(gr)$height <- 25
  graph::nodeRenderInfo(gr)$lWidth <- 37.5
  graph::nodeRenderInfo(gr)$rWidth <- 37.5
  graph::edgeRenderInfo(gr)$weight <- 20

  if (render == TRUE) {
    Rgraphviz::renderGraph(gr)
  }

  invisible(gr)
}
