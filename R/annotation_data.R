#' Map of PTK PamChip IDs to PTK
#'
#' A dataset that contains a map of all 196 PTK PamChip Peptides
#' and the PT Kinases that could phosphorylate them in a human sample
#'
#' @format A dataframe with 991 rows and 6 columns
#' \describe{
#' \item{ID}{The ID of the Peptide on the PTK PamChip}
#' \item{Group}{The Group categorization for the associated kinase from RegPhos 2.0 Annoatation}
#' \item{Family}{The Group categorization for the associated kinase from RegPhos 2.0 Annoatation}
#' \item{Subfamily}{The Subfamily categorization for the associated kinase from RegPhos 2.0 Annoatation. Can be NA}
#' \item{Kinase}{The Kinase tha can potentially phosphorylate this peptide}
#' \item{Gene_Symbol}{The HGNC Gene Symbol for the Kinase}
#' }
"ptk_annotation"

#' Protein-Protein Interactome of humans
#'
#' A complete known map of Protein-Protein interactions in the human
#' interactome
#'
#' @format A dataframe with 137889 rows and 2 columns
#'
#' \describe{
#' \item{from}{The Affecting protein}
#' \item{to}{The Affected Protein}
#' }
#'
"ptk_interactome"

#' Map of STK PamChip IDs to PTK
#'
#' A dataset that contains a map of all 142 PTK PamChip Peptides
#' and the PT Kinases that could phosphorylate them in a human sample
#'
#' @format A dataframe with 991 rows and 6 columns
#' \describe{
#' \item{ID}{The ID of the Peptide on the PTK PamChip}
#' \item{Group}{The Group categorization for the associated kinase from RegPhos 2.0 Annoatation}
#' \item{Family}{The Group categorization for the associated kinase from RegPhos 2.0 Annoatation}
#' \item{Subfamily}{The Subfamily categorization for the associated kinase from RegPhos 2.0 Annoatation. Can be NA}
#' \item{Kinase}{The Kinase tha can potentially phosphorylate this peptide}
#' \item{Gene_Symbol}{The HGNC Gene Symbol for the Kinase}
#' }
"stk_annotation"

#' Protein-Protein Interactome of humans
#'
#' A complete known map of Protein-Protein interactions in the human
#' interactome
#'
#' @format A dataframe with 137889 rows and 2 columns
#'
#' \describe{
#' \item{from}{The Affecting protein}
#' \item{to}{The Affected Protein}
#' }
#'
"stk_interactome"
