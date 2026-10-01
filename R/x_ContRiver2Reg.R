#' @title Continental surface water to Regional [s-1]
#' @name x_ContRiver2Reg
#' @param ScaleName Name of the relevant scale
#' @param SubCompartName Name of the relevant sub-compartment
#' @param RunoffFlow RunoffFlow flow from soil to river [m3.s-1]
#' @param RainOnFreshwater Water flow of rain directly on lake/river [m3.s-1]
#' @param dischargeFRAC Fraction discharge of freshwater between regional and continental scales and vice versa [-]
#' @return River Discharge for scale Continental [s-1]
#' @export
x_ContRiver2Reg <- function(ScaleName, SubCompartName, RunoffFlow, RainOnFreshwater,dischargeFRAC, parent) {
  
  # Calculate a flow from
  # Input often contains data with different keys that are needed for a value
  # Runoff needs to be summed across soils compartment to represent input for river compartment
  data <- ScaleName |>
    tidyr::expand_grid(SubCompartName) |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::filter(SubCompartName  == 'river') |>
    dplyr::full_join(RunoffFlow, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(RainOnFreshwater, by=c("Scale", "SubCompart")) |>
    dplyr::full_join(dischargeFRAC, by="Scale") |>
    dplyr::filter(Scale == 'Continental') |>
    dplyr::summarise(
      Scale = 'Continental',
      SubCompart = 'river',
      Runoff  = sum(Runoff , na.rm=TRUE),
      RainOnFreshwater = sum(RainOnFreshwater, na.rm=TRUE),
      dischargeFRAC = dplyr::first(dischargeFRAC)
    ) |>
    dplyr::mutate(
      x_ContRiver2Reg = (Runoff + RainOnFreshwater) * dischargeFRAC
    ) |>
    dplyr::select(Scale, SubCompart, x_ContRiver2Reg)
  
  # The above calculated from joined with the scaffolding, leftjoin should preserve scaffold
  out <- parent$FromDataAndTo("x_ContRiver2Reg") |>
    dplyr::left_join(data, by=c("from.Scale" = "Scale", "SubCompart" = "SubCompart"))
  
  return(data.frame(out)) 
  
  # switch(ScaleName,
  #   "Continental" = {
  #     switch(SubCompartName,
  #       "river" = {
  #         SumRainRunoff <- sum(all.RunoffFlow$RunoffFlow[all.RunoffFlow$Scale == ScaleName]) +
  #           sum(all.RainOnFreshwater$RainOnFreshwater[all.RainOnFreshwater$Scale == ScaleName])
  #         # River2sea  <- RainOnFreshwater + SumRunoff * (1-dischargeFRAC)
  #         # Lake2River <- LakeFracRiver * River2sea
  # 
  #         return((SumRainRunoff) * dischargeFRAC)
  #       },
  #       return(NA)
  #     )
  #   },
  #   return(NA)
  # )
}
