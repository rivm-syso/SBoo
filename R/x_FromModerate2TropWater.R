#' @title Flow from Moderate to Tropic deep ocean (w3)
#
#' Moderate scale waters flow into Tropical scale through deep ocean currents 
#' based on a set ocean current of r World$fetchData("OceanCurrent")/(60*60*24) m3/day.
#' 
#' @param OceanCurrent Ocean Current between global zones [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Flow [m3.s-1]
#' @export
#' 
x_FromModerate2TropWater <- function ( OceanCurrent, SubCompartName, ScaleName, parent) {
  
  out <- parent$FromDataAndTo("x_FromModerate2TropWater") |>
    dplyr::mutate(
      x_FromModerate2TropWater = dplyr::case_when(
        SubCompart == "sea" ~ 0,
        SubCompart == "deepocean" ~ OceanCurrent, #fix, data in code
        TRUE ~ NA
      )
    )
  return(data.frame(out))
  
           # switch (SubCompartName,
           #         "sea" = {return(0)},
           #         "deepocean" = {OceanCurrent}, #fix, data in code
           #         NA
  # )
  
}
