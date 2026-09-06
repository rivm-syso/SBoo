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
    expand_grid(SubCompartName, SpeciesName) |>
    full_join(Area, by=c("Scale", "SubCompart")) |>
    full_join(Volume, by=c("Scale", "SubCompart")) |>
    full_join(WINDspeed, by=c("Scale")) |>
    filter(SubCompartName %in% c("air")) |>
    parent$states$clipStates() |>
    mutate(
      TAU = dplyr::case_when(
        SubCompart %in% c("air") ~ f_TAU(Area, WINDspeed),
        TRUE ~ NA_real_
      ),
      AirFlow = Volume / TAU
    ) |>
    filter(!is.na(AirFlow)) |>
    select(Scale, SubCompart, AirFlow) |>
    arrange(Scale, SubCompart)
    
  return(data.frame(out))
    
  # if (SubCompartName %in%  c("air")) { #, "cloudwater" should also?!
  #   TAU <- f_TAU(Area, WINDspeed) #Residence time
  #   Volume / TAU
  # } else {
  #   return(NA) #not compartment "air"; not a valid airflow
  # }
}
