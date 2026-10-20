#' @title Continental surface water to Regional [s-1]
#' @name x_ContRiver2Reg
#' @param ScaleName Name of the relevant scale
#' @param SubCompartName Name of the relevant sub-compartment
#' @param RunoffFlow RunoffFlow flow from soil to river [m3.s-1]
#' @param RainOnFreshwater Water flow of rain directly on lake/river [m3.s-1]
#' @param dischargeFRAC Fraction discharge of freshwater between regional and continental scales and vice versa [-]
#' @return River Discharge for scale Continental [s-1]
#' @export
x_ContRiver2Reg <- function(RunoffFlow, RainOnFreshwater, dischargeFRAC, parent) {
  
  out <- parent$FromDataAndTo("x_ContRiver2Reg")

  Runoff <- RunoffFlow |>
    dplyr::group_by(Scale) |>
    dplyr::summarise(allsoils_runoff = sum(RunoffFlow, na.rm = TRUE))
  
  # Runoff needs to be summed across soils compartment to represent input for river compartment
  out <- out |> 
    dplyr::left_join(Runoff, by = c("from.Scale" = "Scale")) |>
    dplyr::left_join(dischargeFRAC, by = c("from.Scale" = "Scale")) |>
    dplyr::left_join(RainOnFreshwater, by = c("from.Scale" = "Scale", "SubCompart")) |>
    
    dplyr::mutate(
      flow = (allsoils_runoff + RainOnFreshwater) * dischargeFRAC
    )

  
  out <- dplyr::select(out, -allsoils_runoff, -dischargeFRAC, -RainOnFreshwater)
  
  return(out) 
  
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
