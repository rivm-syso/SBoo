#' @title Fragmentation of macroplastics to microplastics
#' @name k_MacroFragmentation
#' @description To Be Completed
#' @param kMfrag Falling appart or fragmenting rate constant [s-1]
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_MacroFragmentation [s-1]
#' @export


k_MacroFragmentation <- function (kMfrag,
                                  SubCompartName,
                                  ScaleName,
                                  Regional_and_Continental_deepocean) {
  if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) &
       (
         SubCompartName == "freshwatersediment" |
         SubCompartName == "lakesediment" |
         SubCompartName == "lake" |
         SubCompartName == "river" |
         SubCompartName == "agriculturalsoil" |
         SubCompartName == "othersoil"
       )
  )) {
    return(NA)
  } else if (!is.null(ScaleName) &&
             ScaleName %in% c("Regional", "Continental") &&
             (SubCompartName == "deepocean") &&
             (
               is.na(Regional_and_Continental_deepocean) ||
               isFALSE(Regional_and_Continental_deepocean) ||
               Regional_and_Continental_deepocean == "FALSE"
             )) {
    return(NA)
  } else {
    return(kMfrag)
  }
  
}