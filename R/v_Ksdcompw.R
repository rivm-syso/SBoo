#'@title Dimensionless sed/water PARTITION COEFFICIENT for molecular species specific to sediment compartment
#'@name Ksdcompw
#'@description Dimensionless sediment/water partitioning coefficient dependent on the compartment
#'@param FRACw fraction of water in compartment [-]
#'@param FRACs fraction of soil in compartment [-]
#'@param Kp general sediment water partitioning coefficient [-]
#'@param rhoMatrix density of the matrix [kg m-3]
#'@param Matrix type of compartment considered
#'@return Ksdcompw
#'@export
Ksdcompw <- function(FRACw, FRACs, Kp, rhoMatrix, Matrix, parent, ScaleName){
  
  RHOsolid <- rhoMatrix |> dplyr::filter(Matrix == 'soil') |> dplyr::pull(rhoMatrix)
  out <- ScaleName |>
    tidyr::expand_grid(Matrix) |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::full_join(FRACw, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(FRACs, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(Kp, by="SubCompart") |>
    dplyr::mutate(
      Ksdcompw = dplyr::case_when(
        Matrix == 'sediment' ~ FRACw+FRACs*Kp*RHOsolid/1000,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(Ksdcompw)) |>
    dplyr::select(Scale, SubCompart, Ksdcompw) |>
    dplyr::arrange(Scale, SubCompart)
  
  return(data.frame(out))
  
  # if (Matrix == "sediment") {
  #   RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  #   return(FRACw+FRACs*Kp*RHOsolid/1000) # we need to take care that RHOsolid here can be specific to the compartment compared to generic one used in Kp in relation to Ksw!
  # } else
  #   return(NA)
}
