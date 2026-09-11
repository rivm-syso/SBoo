#'@title Fraction of chemical truly dissolved in water phase fresh water (relevant to Molecular species)
#'@name Frinw
#'@description Everything that is not attached
#'@param Kp soil/water partitioning coefficient [-]
#'@param SUSP Concentration suspended matter in water [mg.L-1]
#'@param COL Concentration of colloidal organic matter in water [mg.L-1]
#'@param Matrix Type of compartment matrix (soil, water, sediment or air)
#'@param FRACa Fraction air in compartment [-]
#'@param FRACw Fraction water in compartment [-]
#'@param FRACs Fraction solids in compartment [-]
#'@param Kacompw Dimensionless air/water partitioning coefficient of original species at compartment temperature [-]
#'@param FRorig_spw Fraction original species in porewater of soil [-]
#'@param rhoMatrix density of the matrix [kg.m-3]
#'@export
FRinw <- function(FRorig_spw, FRACw, FRACa, FRACs, Kp, rhoMatrix, KpCOL, 
                  Kacompw, SUSP, COL, Matrix, parent, ScaleName){
  
  RHOsolid <- rhoMatrix |> dplyr::filter(Matrix == 'soil') |> dplyr::pull(rhoMatrix)
  out <- ScaleName |>
    tidyr::expand_grid(Matrix) |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::full_join(FRACw, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(FRACa, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(FRACs, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(FRorig_spw, by='SubCompart') |>
    dplyr::full_join(Kp, by="SubCompart") |>
    dplyr::full_join(KpCOL, by="SubCompart") |>
    dplyr::full_join(Kacompw, by="Scale") |>
    dplyr::full_join(SUSP, by ="SubCompart") |>
    dplyr::full_join(COL, by="SubCompart") |>
    dplyr::mutate(
      Frinw = dplyr::case_when(
        Matrix == 'water' ~ 1/(1+Kp*SUSP/1000+KpCOL*COL/1000),
        Matrix == 'soil' ~ FRACw/(FRACa*(Kacompw*FRorig_spw)+FRACw+FRACs*Kp*RHOsolid/1000),
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(Frinw)) |>
    dplyr::select(Scale, SubCompart, Frinw) |>
    dplyr::arrange(Scale, SubCompart)
    
  return(data.frame(out))
   
  # RHOsolid <- all.rhoMatrix$rhoMatrix[all.rhoMatrix$SubCompart == "naturalsoil"]
  # switch(Matrix,
  #        "water" = 1/(1+Kp*SUSP/1000+KpCOL*COL/1000),
  #        "soil" = # for soil pore water
  #          FRACw/(FRACa*(Kacompw*FRorig_spw)+FRACw+FRACs*Kp*RHOsolid/1000)
  # )
}
