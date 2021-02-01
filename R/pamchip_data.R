#' @import methods
NULL

#' The Pamchip Data Superclass
#'
#' @slot chip_type character. A string
#'
#' @return an object of class PamchipData
#'
#' @examples
#' TRUE
setClass("PamchipData",
         slots = c(
           ChipType = "character"
         ),
         contains = "VIRTUAL")

setValidity("PamchipData",
            function(object) {
              if(!object@ChipType %in% c("PTK", "STK")) {
                "Invalid @chip_type"
              } else {
                TRUE
              }
            })
