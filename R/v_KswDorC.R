#' @title Ksw
#' @name KswDorC
#' @param Kow
#' @param pKa
#' @param CorgStandard
#' @param ChemClass
#' @param a
#' @param b
#' @param RHOsolid
#' @param Ksw
#' @export
KswDorC <- function (Kow, pKa, CorgStandard, ChemClass, a, b, all.rhoMatrix, Ksw){
  RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  
  if (is.na(Ksw) || Ksw == "NA") { 
    f_Ksw (Kow, pKa, CorgStandard , a, b, ChemClass, RHOsolid, FALSE, Ksw)
  } else return(Ksw)
  
}
