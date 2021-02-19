#' Generate a reduced kianse graph
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
#' @importFrom tibble as_tibble
#'
#' @examples
#' TRUE
render_reduced_kinased_graph <- function(analysis_result, title, render = FALSE) {
  mapped_kinases <- analysis_result$mapped_kinases %>%
    dplyr::select(.data$peptide, .data$kinase) %>%
    dplyr::group_by(.data$peptide) %>%
    dplyr::slice_head(n = 1) %>%
    tibble::deframe()

  net_arcs <- bnlearn::arcs(analysis_result$net)

  mapped_net_arcs <- net_arcs %>%
    tibble::as_tibble() %>%
    mutate(
      from = mapped_kinases[.data$from],
      to = mapped_kinases[.data$to]
    ) %>%
    unique() %>%
    dplyr::filter(.data$from != .data$to) %>%
    as.matrix()

  bn <- bnlearn::empty.graph(unique(c(mapped_net_arcs)))
  bnlearn::arcs(bn, check.cycles = FALSE, check.illegal = FALSE) <- mapped_net_arcs

  network_graph <- bnlearn::graphviz.plot(bn,
                                          shape = "ellipse", layout = "fdp",
                                          main = title)


  gr <- Rgraphviz::layoutGraph(network_graph, attrs = list(graph = list(rankdir = "TB")))
  graph::graphRenderInfo(gr)$fontsize <- 24
  graph::nodeRenderInfo(gr)$fontsize <- 14
  graph::nodeRenderInfo(gr)$fixedsize <- TRUE
  graph::nodeRenderInfo(gr)$height <- 25
  graph::nodeRenderInfo(gr)$lWidth <- 37.5
  graph::nodeRenderInfo(gr)$rWidth <- 37.5
  graph::edgeRenderInfo(gr)$weight <- 20

  if (render == TRUE) {
    Rgraphviz::renderGraph(gr)
  }

  out <- list(bnnet = bn,
              graphviz_plot = gr)

  invisible(out)
}
