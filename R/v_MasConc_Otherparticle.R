# POTENTIALLY OBSOLETE

#' @title MasConc_Otherparticle
#' @name MasConc_Otherparticle
#' @param mConcSusp concentration of ... in 
#' @param mConcCol concentration of ... in 
#' @param SpeciesName yes, its name
#' @return MasConc_Otherparticle, mass concentration of "other" particle before 
#' @export
MasConc_Otherparticle <- function (mConcSusp, mConcCol, SpeciesName){
  dplyr::case_when(
    SpeciesName == "Small" ~ mConcCol,
    SpeciesName == "Large" ~ mConcSusp,
    TRUE ~ NA_real_  #all but Small or Large
  )
}
