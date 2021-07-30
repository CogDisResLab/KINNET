#' @include pamchip_data_ptk.R
#' @import methods stringr
NULL

setMethod(
  "show",
  signature = "PamchipData-PTK",
  definition = function(object) {
    print(
      stringr::str_glue(
        "{object@ChipType} PamChip Data
                        Processed with BioNavigator {object@BioNavigatorVersion}
                        Processing done on {object@ImageAnalysisDate}
                        {nrow(exp_data(object))} Observations in {length(unique(pheno_data(object)$class))} distinct classes
                        "
      )
    )
    invisible(NULL)
  }
)

setGeneric("processed_data", function(object, ...) {
  standardGeneric("processed_data")
})

setMethod("processed_data", "PamchipData-PTK",
          function(object) {
            object@ProcessedData
          })

setGeneric("processed_data<-", function(object, value) {
  standardGeneric("processed_data<-")
})
setMethod("processed_data<-", "PamchipData-PTK",
          function(object, value) {
            object@ProcessedData <- value
            if (validObject(object)) {
              return(object)
            }
          })

setGeneric("pheno_data", function(object, ...) {
  standardGeneric("pheno_data")
})

setMethod("pheno_data", "PamchipData-PTK",
          function(object) {
            object@SampleCharacteristics
          })

setGeneric("exp_data", function(object, ...) {
  standardGeneric("exp_data")
})

setMethod("exp_data", "PamchipData-PTK",
          function(object) {
            object@SampleData
          })

setGeneric("classes", function(object, ...) {
  standardGeneric("classes")
})

setMethod("classes", "PamchipData-PTK",
          function(object) {
            unique(object@SampleCharacteristics$class)
          })


setGeneric("peptides", function(object, ...) {
  standardGeneric("peptides")
})

setMethod("peptides", "PamchipData-PTK",
          function(object) {
            object@PeptideIDs
          })
