#' @title Input or calculated organic carbon partitioning coefficient for the original form..
#' @name KocDorC
#' @description Function that either calculates the Organic Carbon partitioning coefficient from established QSARs based on Kow and pKa
#' Or uses input data for Koc.
#' @param Kow Octanol water partitioning coefficient in data
#' @param a see QSAR table 
#' @param b see QSAR table 
#' @param Koc Organic Carbon partitioning coefficient in
#' @param pKa Dissociation constant of (conjugated) acid (default = 7)
#' @export
KocDorC <- function (Kow, a, b, Koc, ChemClass){
    
  if (is.na(Koc) || Koc == "NA") {
    if (ChemClass == 'particle') {
      return(NA)
    }
    out <- a |>
      dplyr::full_join(b, by="QSAR.ChemClass") |>
      dplyr::filter(QSAR.ChemClass %in% ChemClass) |>
      dplyr::mutate(
        KocDorC = dplyr::case_when(
          QSAR.ChemClass == 'acid' ~ 10^(0.54*log10(Kow)+1.11),
          QSAR.ChemClass == 'base' ~ 10^(0.37*log10(Kow)+1.7),
          QSAR.ChemClass == 'metal' ~ NA,
          QSAR.ChemClass == 'particle' ~ NA,
          TRUE ~ a * Kow^b
        )
      ) |>
      dplyr::pull(KocDorC)
    
    return(out)
  } else return(Koc)
    
  #   if (is.na(Koc) || Koc == "NA") { 
  #     switch(ChemClass,
  #            "acid" = 10^(0.54*log10(Kow)+1.11) ,
  #            "base" = 10^(0.37*log10(Kow)+1.7) ,
  #            "metal" = NA,
  #            "particle" = NA,
  #            #else
  #            {a * Kow^b})
  # 
  # } else return(Koc)
  # 
}