#' Return Parents/Children with a Usable Default
#'
#' @param x A bnlearn network
#' @param n A node
#'
#' @return Parents or children, given the function. NA if no parents or no children
#' @import bnlearn
#'
#' @keywords internal
#' @examples
#' TRUE
default_parent <- function(x, n) {
  p <- bnlearn::parents(x, n)
  if (length(p) == 0) {
    NA
  } else {
    p
  }
}

#' @describeIn default_parent Return a default if children don't exist
default_children <- function(x, n) {
  p <- bnlearn::children(x, n)
  if (length(p) == 0) {
    NA
  } else {
    p
  }
}

#' Generate a Triad of Peptides Given A Network
#'
#' @param arcs a dataframe of arcs between peptides
#' @param net a bnlearn network
#'
#' @return a dataframe with all combinations of parents and children of a peptide as kinases.
#' @import purrr
#'
#' @keywords internal
#'
#' @examples
#' TRUE
generate_triad_links <- function(arcs, net) {
  triad_links <- arcs %>% purrr::map_dfr(~ expand_grid(
    parent = default_parent(net, .x),
    node = .x,
    child = default_children(net, .x)
  )) %>%
    unique

  triad_links
}

#' Filter list of kinases based on the parents and children
#'
#' @param parent Parent node for the given triad
#' @param node The node under consideration
#' @param child Child node for the given triad
#' @param kinase_annotation A dataframe with the peptide to kinase map
#' @param interactome A dataframe with the kinase-kinase interactions
#'
#' @return A dataframe of filtered kinases
#' @import dplyr
#'
#' @keywords internal
#'
#' @examples
#' TRUE
process_parent_child_link <- function(parent, node, child, kinase_annotation, interactome) {
  if (!is.na(parent) & !is.na(child)) {
    parent_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == parent) %>%
      dplyr::pull(.data$Gene_Symbol)

    node_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == node) %>%
      dplyr::pull(.data$Gene_Symbol)

    child_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == child) %>%
      dplyr::pull(.data$Gene_Symbol)

    expanded_pn <- expand.grid(
      from = parent_linked,
      to = node_linked
    ) %>% tibble::as_tibble()

    expanded_nc <- expand.grid(
      from = node_linked,
      to = child_linked
    ) %>% tibble::as_tibble()

    filtered_pn <- dplyr::intersect(expanded_pn, interactome) %>%
      dplyr::rename(parent = .data$from,
                    node = .data$to)
    filtered_nc <- dplyr::intersect(expanded_nc, interactome) %>%
      dplyr::rename(node = .data$from,
                    child = .data$to)

    filtered_links <- dplyr::inner_join(filtered_pn, filtered_nc, by = "node")
    filtered_links
  } else if (is.na(parent) & !is.na(child)) {
    parent_linked <- NA

    node_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == node) %>%
      dplyr::pull(.data$Gene_Symbol)

    child_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == child) %>%
      dplyr::pull(.data$Gene_Symbol)

    expanded_nc <- expand.grid(
      from = node_linked,
      to = child_linked
    ) %>% tibble::as_tibble()

    filtered_nc <- dplyr::intersect(expanded_nc, interactome) %>%
      dplyr::rename(node = .data$from,
                    child = .data$to) %>%
      dplyr::mutate(parent = NA) %>%
      dplyr::select(.data$parent, .data$node, .data$child)

    filtered_links <- filtered_nc
    filtered_links
  } else {
    parent_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == parent) %>%
      dplyr::pull(.data$Gene_Symbol)

    node_linked <- kinase_annotation %>%
      dplyr::filter(.data$ID == node) %>%
      dplyr::pull(.data$Gene_Symbol)

    child_linked <- NA

    expanded_pn <- expand.grid(
      from = parent_linked,
      to = node_linked
    ) %>% tibble::as_tibble()

    filtered_pn <- dplyr::intersect(expanded_pn, interactome) %>%
      dplyr::rename(parent = .data$from,
                    node = .data$to) %>%
      dplyr::mutate(child = NA) %>%
      dplyr::select(.data$parent, .data$node, .data$child)

    filtered_links <- filtered_pn
    filtered_links
  }

}

#' Generate Intersections between kinases
#'
#' @param row row peptide
#' @param column column peptide
#' @param assigned_kinases Kinases assigned
#'
#' @return common kinases beween the two groups
#'
#' @examples
#' TRUE
get_intersections <- function(row, column, assigned_kinases) {
  out <- assigned_kinases[[as.numeric(row)]][[column]]
  out
}

#' Generate Candidate Kinases Based On Links
#'
#' @param peptide Peptide ID
#' @param arcs The arcs in the network
#' @param assigned_kinases The candidate kinases
#'
#' @return common kinases
#' @import purrr
#'
#' @examples
#' TRUE
candidate_kinases <- function(peptide, arcs, assigned_kinases) {
  arcs %>%
    rownames_to_column("row") %>%
    filter(peptide == .data$parent | peptide == .data$node | peptide == .data$child) %>%
    mutate(column = case_when(
      peptide == .data$parent ~ "parent",
      peptide == .data$node ~ "node",
      peptide == .data$child ~ "child"
    )) %>%
    select(.data$row, .data$column) %>%
    purrr::pmap(~ get_intersections(..1, ..2, assigned_kinases)) %>%
    purrr::reduce(intersect)
}

#' Assign kinases, given a network and a chip type
#'
#' @param network A network output from bnlearn
#' @param chiptype Either PTK or STK
#'
#' @return A dataframe with upstream kinases assigned to the peptide
#' @export
#'
#' @import dplyr bnlearn purrr tidyr
#'
#' @examples
#' TRUE
assign_kinases <- function(network, chiptype = "PTK") {

  # Assign the appropriate annotation
  if (!chiptype %in% c("PTK", "STK")) {
    stop("Incorrect Chip Type specified")
  }

  if (chiptype == "PTK") {
    annotation <- JustinKinomeModelling::ptk_annotation
    interactome <- JustinKinomeModelling::ptk_interactome
  } else if (chiptype == "STK") {
    annotation <- JustinKinomeModelling::stk_annotation
    interactome <- JustinKinomeModelling::stk_interactome
  }






  # Load and filter appropriate data
  subset_kinase_annotation <- annotation %>%
    dplyr::filter(.data$ID %in% unique(bnlearn::arcs(network)))

  subset_interactome <- interactome %>%
    dplyr::filter(.data$from %in% subset_kinase_annotation$Gene_Symbol |
                    .data$to %in% subset_kinase_annotation$Gene_Symbol)



  # Get the arcs and the unique peptides
  network_arcs <- bnlearn::arcs(network)
  network_peptides <- unique(c(network_arcs))
  names(network_peptides) <- network_peptides

  triad_links <- generate_triad_links(network_arcs, network)

  # Generate a list of all possible kinases
  assigned_links <- purrr::pmap(triad_links, process_parent_child_link,
                                kinase_annotation = subset_kinase_annotation,
                                interactome = subset_interactome)

  # Assign Kinases
  kinase_assignment <- purrr::map(network_peptides, ~ candidate_kinases(.x, triad_links, assigned_links))

  # Generate a DF
  kinase_assignment_df <-
    purrr::map_dfr(names(kinase_assignment), function(x)
      tidyr::expand_grid(peptide = x, kinase  = kinase_assignment[[x]])) %>%
    dplyr::rename(ID = .data$peptide,
                  Gene_Symbol = .data$kinase)

  updated_probabilities <- update_probability_matrix(chiptype, kinase_assignment_df)

  mapped_kinases <- updated_probabilities %>%
    dplyr::group_by(.data$peptide) %>%
    dplyr::filter(.data$logged.foldchange == max(.data$logged.foldchange))

  out <- list(
    net = network,
    probability_df = updated_probabilities,
    kinase_assignment = kinase_assignment_df,
    mapped_kinases = mapped_kinases
  )

  out
}
