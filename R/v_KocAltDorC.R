#' @title Input or calculated organic carbon partitioning coefficient for the alternative form.
#' @name KocAltDorC
#' @description Function that either calculates the Organic Carbon partitioning coefficient from established QSARs based on Kow and pKa
#' Or uses input data for KocAlt.
#' @param Kow Octanol water partitioning coefficient in data
#' @param a see QSAR table 
#' @param b see QSAR table 
#' @param KocAlt Organic Carbon partitioning coefficient in data
#' @param pKa Dissociation constant of (conjugated) acid (default = 7)

#' @export
KocAltDorC <- function (Kow, a, b, pKa, KocAlt, ChemClass){

    if (is.na(KocAlt) || KocAlt == "NA") {
      if (is.na(pKa) || pKa == "NA"){
        pKa <- 7
        warning("KswDorC: pKa is needed but missing, setting pKa=7", call. = FALSE)
      }
      
      out <- a |>
        dplyr::full_join(b, by="QSAR.ChemClass") |>
        dplyr::filter(QSAR.ChemClass %in% ChemClass) |>
        dplyr::mutate(
          KocAltDorC = dplyr::case_when(
            QSAR.ChemClass == 'acid' ~ 10^(0.11*log10(Kow)+1.54),
            QSAR.ChemClass == 'base' ~ 10^(pKa^0.65*(Kow/(1+Kow))^0.14),
            TRUE ~ a * Kow^b
            )
          ) |>
        dplyr::pull(KocAltDorC)
        
        return(out)
      # switch(ChemClass,
      #        "acid" = 10^(0.11*log10(Kow)+1.54) ,
      #        "base" = 10^(pKa^0.65*(Kow/(1+Kow))^0.14),
      #        #else
      #        {a * Kow^b}      
      
    
  } else return(KocAlt)
  
}
