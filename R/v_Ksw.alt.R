#' @title Ksw.alt
#' @name Ksw.alt
#' @param Kow
#' @param pKa
#' @param CorgStandard
#' @param ChemClass
#' @param a
#' @param b
#' @param rhoMatrix Density of the matrix used for getting RHOsolid [kg.m-3]
#' @param KswDorC
#' @export
Ksw.alt <- function (Kow, pKa, CorgStandard, ChemClass, a, b, all.rhoMatrix, KswDorC){
  RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  f_Ksw(Kow=Kow, 
        pKa=pKa, 
        CorgStandard=CorgStandard , 
        a=a, 
        b=b, 
        ChemClass=ChemClass, 
        RHOsolid=RHOsolid, 
        alt_form=TRUE, 
        Ksw_orig=KswDorC)
}
