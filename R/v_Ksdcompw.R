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
Ksdcompw <- function(FRACw, FRACs, Kp, rhoMatrix, Matrix, parent, SpeciesName, ScaleName){
  
  RHOsolid <- rhoMatrix |> filter(Matrix == 'soil') |> pull(rhoMatrix)
  out <- ScaleName |>
    expand_grid(Matrix, SpeciesName) |>
    parent$states$clipStates() |>
    full_join(FRACw, by=c("Scale", "SubCompart")) |>
    full_join(FRACs, by=c("Scale", "SubCompart")) |>
    full_join(Kp, by="SubCompart") |>
    mutate(
      Ksdcompw = dplyr::case_when(
        Matrix == 'sediment' ~ FRACw+FRACs*Kp*RHOsolid/1000,
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(Ksdcompw)) |>
    select(Scale, SubCompart, Ksdcompw) |>
    arrange(Scale, SubCompart)
    
  return(data.frame(out))
  
  # if (Matrix == "sediment") {
  #   RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  #   return(FRACw+FRACs*Kp*RHOsolid/1000) # we need to take care that RHOsolid here can be specific to the compartment compared to generic one used in Kp in relation to Ksw!
  # } else
  #   return(NA)
}
