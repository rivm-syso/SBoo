#' @title Dimensionless soil/water Partition Coefficient of molecular species specific to soil subcompartment
#' @name Kscompw
#' @description Dimensionless soil/water partition coefficient specific to subcompartment
#' @param FRACw fraction of water in compartment [-]
#' @param FRACa fraction of air in compartment [-]
#' @param Kp subcompartment/water partitioning coefficient [-]
#' @param FRorig_spw fraction of original species in soil pore water [-]
#' @param rhoMatrix density of the matrix [-]
#' @param Matrix type of compartment
#' @return Kscompw
#' @export
Kscompw <- function(FRACw, FRACa, Kacompw, FRorig_spw, Kp, rhoMatrix, Matrix, parent, SpeciesName, ScaleName) {
  
  RHOsolid <- rhoMatrix |> filter(Matrix == 'soil') |> pull(rhoMatrix)
  out <- ScaleName |>
    expand_grid(Matrix, SpeciesName) |>
    parent$states$clipStates() |>
    full_join(FRACw, by=c("Scale", "SubCompart")) |>
    full_join(FRACa, by=c("Scale", "SubCompart")) |>
    full_join(Kp, by="SubCompart") |>
    full_join(FRorig_spw, by = "SubCompart") |>
    full_join(Kacompw, by = "Scale") |>
    mutate(
      Kscompw = dplyr::case_when(
        Matrix == 'soil' ~ FRACa * (Kacompw * FRorig_spw) + FRACw + (1 - FRACa - FRACw) * Kp * RHOsolid / 1000,
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(Kscompw)) |>
    select(Scale, SubCompart, Kscompw) |>
    arrange(Scale, SubCompart)
  
  return(data.frame(out))
  
  # if (Matrix == "soil") {
  #   RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  #   return(FRACa * (Kacompw * FRorig_spw) + FRACw + (1 - FRACa - FRACw) * Kp * RHOsolid / 1000)
  #   # we need to take care that RHOsolid above can be specific to the compartment
  # } else {
  #   return(NA)
  # }
}
