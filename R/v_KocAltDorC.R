#' @title Ksw
#' @name KocDorC
#' @description soil-water partitioning coefficient for organic colloids
#' @param a see QSAR table 
#' @param b see QSAR table 
#' @param rhoMatrix density of the matrix [kg/m3]
#' @param pKa Dissociation constant of (conjugated) acid (default = 70
#' @param Ksw soil water partitioning coefficient in data
#' @export
KocAltDorC <- function (Kow, a, b, KocAlt){

    if (is.na(KocAlt) || Koc == "NA") { 
    return(a * Kow^b)
    
  } else return(KocAlt)
  
}
