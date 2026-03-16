#' @title Ksw
#' @name KocDorC
#' @description soil-water partitioning coefficient for organic colloids
#' @param a see QSAR table 
#' @param b see QSAR table 
#' @param rhoMatrix density of the matrix [kg/m3]
#' @param pKa Dissociation constant of (conjugated) acid (default = 70
#' @param Ksw soil water partitioning coefficient in data
#' @export
KocAltDorC <- function (Kow, a, b, pKa, KocAlt){

    if (is.na(KocAlt) || KocAlt == "NA") { 
      if (is.na(pKa) || pKa == "NA"){
        pKa <- 7
        warning("KswDorC: pKa is needed but missing, setting pKa=7", call. = FALSE)
      }
      switch(ChemClass,
             "acid" = 10^(0.11*log10(Kow)+1.54) ,
             "base" = 10^(pKa^0.65*(Kow/(1+Kow))^0.14),
             #else
             {a * Kow^b}      
      )
    
  } else return(KocAlt)
  
}
