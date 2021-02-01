test_that("PamchipData_PTK Process a file into a usable structure", {
  #loading the test data
  dataset <- "PTK-Test.txt"

  test_object <- PamchipData_PTK(dataset)

  #generate initial reference variables
  input_data <- readLines(dataset, warn = F, encoding = "UTF-8")
  t_id_line <- which(stringr::str_detect(input_data, "ID"))
  t_ref_end <- max(which(stringr::str_detect(input_data, "REF")))
  t_metadata_line <- stringr::str_c(input_data[1], "\n")
  t_characteristic_data <- stringr::str_c(input_data[2:(t_id_line-1)], collapse = "\n")
  t_ref_data <- stringr::str_c(input_data[(t_id_line + 1):t_ref_end], collapse = "\n")
  t_sample_data <- stringr::str_c(input_data[(t_ref_end + 1):length(input_data)], collapse = "\n")

  #read and split the data
  input_data <- readLines(dataset, warn = F, encoding = "UTF-8")
  id_line <- which(stringr::str_detect(input_data, "ID"))
  ref_end <- max(which(stringr::str_detect(input_data, "REF")))
  metadata_line <- stringr::str_c(input_data[1], "\n")
  characteristic_data <- stringr::str_c(input_data[2:(id_line-1)], collapse = "\n")
  ref_data <- stringr::str_c(input_data[(id_line + 1):ref_end], collapse = "\n")
  sample_data <- stringr::str_c(input_data[(ref_end + 1):length(input_data)], collapse = "\n")

  #checkpoint 1
  expect_identical(t_id_line, id_line)
  expect_identical(t_ref_end, ref_end)
  expect_identical(t_metadata_line, metadata_line)
  expect_identical(t_characteristic_data, characteristic_data)
  expect_identical(t_ref_data, ref_data)
  expect_identical(t_sample_data, sample_data)

  #generate refrence metadata
  t_metadata <- readr::read_tsv(t_metadata_line, col_names = F)
  t_BioNavigatorVersion <- t_metadata$X2
  t_ImageAnalysisDate <- t_metadata$X8
  t_QuantitationType <- t_metadata$X5

  #processing metadata
  metadata <- readr::read_tsv(metadata_line, col_names = F)
  BioNavigatorVersion <- metadata$X2
  ImageAnalysisDate <- metadata$X8
  QuantitationType <- metadata$X5

  #checkpoint 2
  expect_identical(t_metadata, metadata)
  expect_identical(t_BioNavigatorVersion, BioNavigatorVersion)
  expect_identical(t_ImageAnalysisDate, ImageAnalysisDate)
  expect_identical(t_QuantitationType, QuantitationType)

  #generate refrence sample characteristics
  t_SampleCharacteristics <- readr::read_tsv(t_characteristic_data, col_names = F) %>%
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

  #process sample characteristics
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

  #generate refrence RefData
  t_RefData <- readr::read_tsv(t_ref_data, col_names = F) %>%
    dplyr::mutate(ID = str_c(gsub("#", "", .data$X1), seq_along(.data$X1), sep = "_")) %>%
    dplyr::select(-.data$X1, -.data$X2) %>%
    dplyr::select(.data$ID, everything()) %>%
    tibble::column_to_rownames("ID") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = str_c("S", stringr::str_pad(seq_along(.data$REF_1), 5, pad = "0"))) %>%
    dplyr::select(.data$sample, everything())


  #process refrence data
  RefData <- readr::read_tsv(ref_data, col_names = F) %>%
    dplyr::mutate(ID = str_c(gsub("#", "", .data$X1), seq_along(.data$X1), sep = "_")) %>%
    dplyr::select(-.data$X1, -.data$X2) %>%
    dplyr::select(.data$ID, everything()) %>%
    tibble::column_to_rownames("ID") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::mutate(sample = str_c("S", stringr::str_pad(seq_along(.data$REF_1), 5, pad = "0"))) %>%
    dplyr::select(.data$sample, everything())

  #generate refrence sample data
  t_SampleData <- readr::read_tsv(t_sample_data, col_names = F) %>%
    dplyr::select(-.data$X2) %>%
    tibble::column_to_rownames("X1") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dummy = "dummy", sample = str_c("S", stringr::str_pad(seq_along(.data$dummy), 5, pad = "0"))) %>%
    dplyr::select(-.data$dummy) %>%
    dplyr::select(.data$sample, everything())

  #process sample data
  SampleData <- readr::read_tsv(sample_data, col_names = F) %>%
    dplyr::select(-.data$X2) %>%
    tibble::column_to_rownames("X1") %>%
    t %>%
    tibble::as_tibble() %>%
    dplyr::mutate(dummy = "dummy", sample = str_c("S", stringr::str_pad(seq_along(.data$dummy), 5, pad = "0"))) %>%
    dplyr::select(-.data$dummy) %>%
    dplyr::select(.data$sample, everything())

  #generate refrence peptide list
  t_PeptideIDs <- colnames(t_SampleData)[2:ncol(t_SampleData)]

  #process peptide list
  PeptideIDs <- colnames(SampleData)[2:ncol(SampleData)]

  #checkpoint 3
  expect_identical(t_SampleCharacteristics, SampleCharacteristics)
  expect_identical(t_RefData, RefData)
  expect_identical(t_SampleData, SampleData)
  expect_identical(t_PeptideIDs, PeptideIDs)


  #create a new object
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
})
