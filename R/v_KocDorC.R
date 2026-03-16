#' @title Ksw
#' @name KocDorC
#' @description soil-water partitioning coefficient for organic colloids
#' @param a see QSAR table 
#' @param b see QSAR table 
#' @param rhoMatrix density of the matrix [kg/m3]
#' @param pKa Dissociation constant of (conjugated) acid (default = 70
#' @param Ksw soil water partitioning coefficient in data
#' @export
KocDorC <- function (Kow, a, b, Koc){

    if (is.na(Koc) || Koc == "NA") { 
      switch(ChemClass,
             "acid" = 10^(0.54*log10(Kow)+1.11) ,
             "base" = 10^(0.37*log10(Kow)+1.7) ,
             "metal" = NA,
             #"particle" = stop("Ksw Should be in the data"),
             "particle" = NA,
             #else
             {a * Kow^b})

  } else return(Koc)
  
}