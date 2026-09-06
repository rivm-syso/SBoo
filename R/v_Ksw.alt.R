#' @title Ksw.alt
#' @name Ksw.alt soil water partitioning coefficient for colloids
#'@param ChemClass Class of chemical, see QSAR table (REACH, 2012)
#'@param Kow octanol water partitioning coefficient [-]
#'@param CorgStandard Standard mass FRACTION organic carbon in soil/sediment [-]
#'@param a see QSAR table 
#'@param b see QSAR table 
#'@param rhoMatrix density of the matrix [kg/m3]
#'@param pKa Dissociation constant of (conjugated) acid (default = 7
#'@param KswDorC soil water partitioning coefficient for colloids [-]
#' @export
Ksw.alt <- function (KocAltDorC, CorgStandard, rhoMatrix){
  out <- rhoMatrix |>
    filter(Matrix == 'soil') |>
    mutate(
      Ksw.alt = KocAltDorC*CorgStandard * rhoMatrix / 1000
    ) |>
    pull(Ksw.alt)
  return(out)
  
  # RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  # return(KocAltDorC*CorgStandard * RHOsolid / 1000)
}
