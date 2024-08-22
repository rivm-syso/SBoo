#' @title ConcentrationModule
#' @description module for calculation of concentrations from masses
#' @import R6
#' @export
ConcentrationModule <- 
  R6::R6Class("ConcentrationModule",
      public = list(
        initialize = function(TheCore, input){
          private$MyCore <- TheCore
          
        }
      ),
      
      private = list(
        MyCore = NULL
      )
  )