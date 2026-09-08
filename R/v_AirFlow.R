#' @title AirFlow
#' @name AirFlow
#' @description Computes the airflow within a compartment
#' @param Volume volume in [m3]
#' @param Area area of the compartment in [m2]
#' @param WINDspeed within the compartment [m.s-1]
#' @param SubCompartName for which the calculation is executed
#' @return AirFlow, from one scale to another 
#' @export
AirFlow <- function (Volume, Area, WINDspeed, SubCompartName, parent, SpeciesName, ScaleName){
  
  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    dplyr::full_join(Area, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(Volume, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(WINDspeed, by=c("Scale")) |>
    dplyr::filter(SubCompartName %in% c("air")) |>
    parent$states$clipStates() |>
    dplyr::mutate(
      TAU = dplyr::case_when(
        SubCompart %in% c("air") ~ f_TAU(Area, WINDspeed),
        TRUE ~ NA_real_
      ),
      AirFlow = Volume / TAU
    ) |>
    dplyr::filter(!is.na(AirFlow)) |>
    dplyr::select(Scale, SubCompart, AirFlow) |>
    dplyr::arrange(Scale, SubCompart)
    
  return(data.frame(out))
    
  # if (SubCompartName %in%  c("air")) { #, "cloudwater" should also?!
  #   TAU <- f_TAU(Area, WINDspeed) #Residence time
  #   Volume / TAU
  # } else {
  #   return(NA) #not compartment "air"; not a valid airflow
  # }
}
