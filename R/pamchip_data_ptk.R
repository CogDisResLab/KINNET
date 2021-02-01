#' @include pamchip_data.R
#' @import methods
NULL


#' A representation of the PamChip PTK Data
#'
#' @slot BioNavigatorVersion character. A string indicating the version of BioNavigator that generated the dataset
#' @slot ImageAnalysisDate character. The date the analysis was conducted on.
#' @slot PamGridVersion character. PamGrid version the chip was run on.
#' @slot QuantitationType character. The kind of qunatitation analysis performed.
#' @slot SampleData tbl_df. A tibble with the observed activity on each peptide
#' @slot SampleCharacteristics tbl_df. A tibble with the characteristics of each sample
#' @slot RefData tbl_df.
#' @slot PeptideIDs character.
#' @slot ProcessedData tbl_df.
#' @slot DataProcessDate character. The date when the data was processed
#'
#' @return An object of class PamchipData-PTK
#' @export
#'
#' @import tibble readr stringr dplyr
#' @importFrom rlang .data
#'
#' @examples
#' TRUE
setClass("PamchipData-PTK",
         slots = c(
           BioNavigatorVersion = "character",
           ImageAnalysisDate = "character",
           PamGridVersion = "character",
           QuantitationType = "character",
           RefData = "tbl_df",
           SampleData = "tbl_df",
           SampleCharacteristics = "tbl_df",
           PeptideIDs = "character",
           ProcessedData = "tbl_df",
           DataProcessDate = "character"
         ),
         contains = "PamchipData",
         prototype = list(
           ChipType = "PTK"
         ))

setValidity("PamchipData-PTK",
            function(object) {
              TRUE
            })


#' Process a file into a usable structure
#'
#' @param dataset character. Path to a file that holds the output of BioNavigator
#'
#' @return An object of class PamchipData-PTK
#' @export
#'
#' @examples
#' TRUE
PamchipData_PTK <- function(dataset) {
  if (!file.exists(dataset)) {
    stop(stringr::str_glue("Dataset {dataset} does not exist."))
  }

  # Read and split data
  input_data <- readLines(dataset, warn = F, encoding = "UTF-8")
  id_line <- which(stringr::str_detect(input_data, "ID"))
  ref_end <- max(which(stringr::str_detect(input_data, "REF")))
  metadata_line <- stringr::str_c(input_data[1], "\n")
  characteristic_data <- stringr::str_c(input_data[2:(id_line-1)], collapse = "\n")
  ref_data <- stringr::str_c(input_data[(id_line + 1):ref_end], collapse = "\n")
  sample_data <- stringr::str_c(input_data[(ref_end + 1):length(input_data)], collapse = "\n")



  # Process Metadata
  metadata <- readr::read_tsv(metadata_line, col_names = F)
  BioNavigatorVersion <- metadata$X2
  ImageAnalysisDate <- metadata$X8
  QuantitationType <- metadata$X5

  # Process Sample Characteristics
  SampleCharacteristics <- readr::read_tsv(characteristic_data, col_names = F) %>%
    dplyr::select(-.data$X1) %>%
    column_to_rownames("X2") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::rename_with(stringr::str_to_lower) %>%
    dplyr::rename_with(function(n) gsub(" ", "_", n)) %>%
    dplyr::rename(class = .data$sample_name,
                  exposure = .data$exposure_time) %>%
    dplyr::mutate(sample = str_c("S", stringr::str_pad(seq_along(.data$barcode), 5, pad = "0"))) %>%
    dplyr::select(.data$sample, everything())

  # Process Reference Data
  RefData <- readr::read_tsv(ref_data, col_names = F) %>%
    dplyr::mutate(ID = str_c(gsub("#", "", .data$X1), seq_along(.data$X1), sep = "_")) %>%
    dplyr::select(-.data$X1, -.data$X2) %>%
    dplyr::select(.data$ID, everything()) %>%
    tibble::column_to_rownames("ID") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = str_c("S", stringr::str_pad(seq_along(.data$REF_1), 5, pad = "0"))) %>%
    dplyr::select(.data$sample, everything())


  # Process Sample Data
  SampleData <- readr::read_tsv(sample_data, col_names = F) %>%
    dplyr::select(-.data$X2) %>%
    tibble::column_to_rownames("X1") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dummy = "dummy",
                  sample = str_c("S", stringr::str_pad(seq_along(.data$dummy), 5, pad = "0"))) %>%
    dplyr::select(-.data$dummy) %>%
    dplyr::select(.data$sample, everything())

  # Get Peptide IDs
  PeptideIDs <- colnames(SampleData)[2:ncol(SampleData)]

  # Create a new object

  chipdata <- new("PamchipData-PTK",
                  BioNavigatorVersion = BioNavigatorVersion,
                  ImageAnalysisDate = ImageAnalysisDate,
                  PamGridVersion = BioNavigatorVersion,
                  QuantitationType = QuantitationType,
                  RefData = RefData,
                  SampleData = SampleData,
                  SampleCharacteristics = SampleCharacteristics,
                  PeptideIDs = PeptideIDs,
                  ProcessedData = tibble(),
                  DataProcessDate = format(Sys.time(), "%a, %d %b %Y %T %Z", tz = "GMT")
  )

  chipdata

}
