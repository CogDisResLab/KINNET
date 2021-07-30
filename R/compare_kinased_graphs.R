#' Add a node safely to a bn object
#'
#' This function takes a bn object and a node and adds it if it's not already in the nodeset.
#'
#' @param net A bn object
#' @param node A node label
#'
#' @return A bn object with the node added
#'
#' @importFrom bnlearn add.node nodes
#'
#' @examples
#' TRUE
safely_add_node <- function(net, node) {
  if (node %in% bnlearn::nodes(net)) {
    net
  } else {
    n <- bnlearn::add.node(net, node)
    n
  }
}

#' Compare two kinased graphhs
#'
#' This functions takes two given graphs and compares them to see the changes from one to the other
#'
#' @param reference The reference network output from assign_kinases()
#' @param comparison The comparison network output from assign_kinases()
#' @param ref_name Name for the reference network
#' @param comp_name Name for the comparison network
#' @param render logical. Whether a graph should be rendered or not.
#'
#' @return A graph object
#' @export
#'
#' @importFrom bnlearn nodes graphviz.compare
#' @importFrom graph nodeRenderInfo edgeRenderInfo graphRenderInfo
#' @importFrom Rgraphviz layoutGraph renderGraph
#' @importFrom stringr str_glue
#'
#' @examples
#' TRUE
compare_kinased_graphs <-
  function(reference,
           comparison,
           ref_name = "Reference",
           comp_name = "Comparison",
           render = FALSE) {
    pair <- equalize_kinase_graphs(reference, comparison)

    title <-
      stringr::str_glue("Comparative Network between {ref_name} and {comp_name}")


    comp_graph <- bnlearn::graphviz.compare(
      pair$reference,
      pair$comparison,
      diff.args = list(
        tp.col = "blue",
        tp.lty = 2,
        fp.col = "darkgreen",
        fp.lty = 1,
        fn.col = "red",
        fn.lty = 1
      ),
      shape = "ellipse",
      layout = "fdp"
    )

    gr <-
      Rgraphviz::layoutGraph(comp_graph[[2]], attrs = list(graph = list(rankdir = "TB", main = title)))
    graph::graphRenderInfo(gr)$fontsize <- 24
    graph::nodeRenderInfo(gr)$fontsize <- 14
    graph::nodeRenderInfo(gr)$fixedsize <- TRUE
    graph::nodeRenderInfo(gr)$height <- 25
    graph::nodeRenderInfo(gr)$lWidth <- 37.5
    graph::nodeRenderInfo(gr)$rWidth <- 37.5
    graph::edgeRenderInfo(gr)$weight <- 20
    graph::edgeRenderInfo(gr)$col <-
      graph::edgeRenderInfo(comp_graph[[2]], "col")

    if (render == TRUE) {
      Rgraphviz::renderGraph(gr)
    }

    invisible(gr)

  }
