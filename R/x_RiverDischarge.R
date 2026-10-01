#' @title RiverDischarge [s-1]
#' @name x_RiverDischarge
#' @param RunoffFlow RunoffFlow flow from soil to river [m3.s-1]
#' @param RainOnFreshwater Water flow of rain directly on lake/river [m3.s-1] 
#' @param dischargeFRAC Fraction discharge of freshwater between regional and continental scales and vice versa [-]
#' @param x_ContRiver2Reg Flow from continental to regional river water [s-1] 
#' @param ScaleName Name of the scale of the box at hand
#' @param SubCompartName Name of the subcompartment of the box at hand
#' @return River Discharge [s-1]
#' @export
x_RiverDischarge <- function (RunoffFlow, RainOnFreshwater, 
                              dischargeFRAC, x_ContRiver2Reg, 
                              ScaleName, SubCompartName, parent){
  
  out <- parent$FromDataAndTo("x_RiverDischarge")
  
  data <- ScaleName |>
    tidyr::expand_grid(SubCompartName) |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::filter(SubCompartName  == 'river') |>
    dplyr::full_join(RunoffFlow, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(RainOnFreshwater, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(dischargeFRAC, by="Scale") |>
    dplyr::filter(Scale %in% out$Scale) |>
    dplyr::group_by(Scale) |>
    
    dplyr::summarise(
      SubCompart = 'river',
      Runoff  = sum(Runoff , na.rm=TRUE),
      RainOnFreshwater = sum(RainOnFreshwater, na.rm=TRUE),
      dischargeFRAC = dplyr::first(dischargeFRAC)
    ) |>
    dplyr::ungroup() |>
    dplyr::mutate(
      x_RiverDischarge = (Runoff + RainOnFreshwater) * (1-dischargeFRAC)
    ) |>
    dplyr::select(Scale, SubCompart, x_RiverDischarge)
  
  out <- dplyr::left_join(out, data, by = c("from.SubCompart" = "SubCompart", "Scale" = "Scale"))
  
  return(data.frame(out))
  
  # x_ContRiver2Reg <- sum(all.x_ContRiver2Reg$flow) #sum to force an atomic number ?
  # SumRainRunoff <- sum(all.RunoffFlow$RunoffFlow[all.RunoffFlow$Scale == ScaleName]) +
  #   sum(all.RainOnFreshwater$RainOnFreshwater[all.RainOnFreshwater$Scale == ScaleName])
  # 
  # switch (SubCompartName,
  #         "river" = {
  #           if(ScaleName == "Continental"){
  #             return((SumRainRunoff) * (1-dischargeFRAC))
  #           } 
  #           if(ScaleName == "Regional"){
  #             return((SumRainRunoff + x_ContRiver2Reg) * (1-dischargeFRAC))
  #           } 
  #           else NA
  #         },
  #         NA
  # )
}
