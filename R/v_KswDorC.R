#' @title Ksw
#' @name KswDorC
#' @description soil-water partitioning coefficient for organic colloids
#' @param ChemClass Class of chemical, see QSAR table (REACH, 2012)
#' @param Kow octanol water partitioning coefficient [-]
#' @param CorgStandard Standard mass FRACTION organic carbon in soil/sediment [-]
#' @param a see QSAR table 
#' @param b see QSAR table 
#' @param rhoMatrix density of the matrix [kg/m3]
#' @param pKa Dissociation constant of (conjugated) acid (default = 70
#' @param Ksw soil water partitioning coefficient in data
#' @export
KswDorC <- function (KocDorC, CorgStandard, rhoMatrix, Ksw, ChemClass){
  
  if (is.na(Ksw) || Ksw == "NA") { 
    out <- rhoMatrix |>
      filter(Matrix == 'soil') |>
      expand_grid(data.frame(ChemClass)) |>
      mutate(
        KswDorC = dplyr::case_when(
          ChemClass == 'metal' ~ NA,
          ChemClass == 'particle' ~ NA,
          TRUE ~ KocDorC*CorgStandard * rhoMatrix / 1000
          )
        ) |>
      pull(KswDorC)
    
    if (out == NA && ChemClass == 'metal') {
      stop("Ksw Should be in the data")
    }
    
    return(out)
  } else return(Ksw)
  
 #  RHOsolid <- rhoMatrix$rhoMatrix[rhoMatrix$SubCompart == "naturalsoil"]
 #  
 #  if (is.na(Ksw) || Ksw == "NA") { 
 #    switch(ChemClass,
 #    #        "acid" = 10^(0.54*log10(Kow)+1.11) ,
 #    #        "base" = 10^(0.37*log10(Kow)+1.7) ,
 #           "metal" = stop("Ksw Should be in the data"),
 #           #"particle" = stop("Ksw Should be in the data"),
 #           "particle" = NA,
 #           #else
 #           { KocDorC*CorgStandard * RHOsolid / 1000})
 #    
 #  }
 #   
 #    
 # else return(Ksw)
  
}
