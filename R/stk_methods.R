#' @include pamchip_data_stk.R
#' @import methods stringr
NULL

setMethod(
  "show",
  signature = "PamchipData-STK",
  definition = function(object) {
    print(
      stringr::str_glue("{object@ChipType} PamChip Data
                        Processed with BioNavigator {object@BioNavigatorVersion}
                        Processing done on {object@ImageAnalysisDate}
                        {nrow(exp_data(object))} Observations in {length(unique(pheno_data(object)$class))} distinct classes
                        ")
    )
    invisible(NULL)
  }
)

#' STK Data Accessor Functions
#'
#' These functions provide a variety of data setters and getter for the Pamchip objects.
#'
#' @details These functions allow you to get and set the slots of the object.
#'
#' \itemize{
#' \item processed_data accesses the processed and transformed data
#' \item pheno_data gives the sample characteristics
#' \item exp_data access the actual intensity values
#' \item classes gives the unique classes in the chip
#' \item peptides gives the reference list of peptides on the chip
#' }
#'
#' @param chipdata an object of class PamchipData-STK or PamchipData-PTK
#' @param ... Currently unused
#'
#' @return The requested object
#'
#' @name stk_accessors
NULL


#' @export
#' @rdname stk_accessors
#'
setGeneric("pheno_data", function(chipdata,...) {
  standardGeneric("pheno_data")
})

#' @export
#' @rdname stk_accessors
#'
setMethod("pheno_data", "PamchipData-STK",
          function(chipdata) {
            chipdata@SampleCharacteristics
          })

#' @export
#' @rdname stk_accessors
#'
setGeneric("exp_data", function(chipdata,...) {
  standardGeneric("exp_data")
})

#' @export
#' @rdname stk_accessors
#'
setMethod("exp_data", "PamchipData-STK",
          function(chipdata) {
            chipdata@SampleData
          })

#' @export
#' @rdname stk_accessors
#'
setGeneric("classes", function(chipdata,...) {
  standardGeneric("classes")
})

#' @export
#' @rdname stk_accessors
#'
setMethod("classes", "PamchipData-STK",
          function(chipdata) {
            unique(chipdata@SampleCharacteristics$class)
          })

#' @export
#' @rdname stk_accessors
#'
setGeneric("peptides", function(chipdata,...) {
  standardGeneric("peptides")
})

#' @export
#' @rdname stk_accessors
#'
setMethod("peptides", "PamchipData-STK",
          function(chipdata) {
            chipdata@PeptideIDs
          })
