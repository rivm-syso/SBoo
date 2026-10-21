#' @title Mixing of upper and deeper sea water layers (w2 and w3)
#' @name x_FromModerate2ContWater
#' @param OceanCurrent [m3.s-1]
#' @param SubCompartName name of the subcompartment of the box at hand
#' @param ScaleName name of the scale of the box at hand
#' @return Advection sea - deepocean
#' @export
#' 
x_FromModerate2ContWater <- function (Volume, TAUsea,
                                     x_RegSea2Cont, parent) {
  
  out <- parent$FromDataAndTo("x_FromModerate2ContWater") |>
    dplyr::left_join(Volume, by=c("to.Scale" = "Scale", "SubCompart" = "SubCompart")) |>
    dplyr::left_join(TAUsea, by=c("to.Scale" = "Scale")) |>
    dplyr::left_join(x_RegSea2Cont |> dplyr::select(to.Scale, flow)) |>
    dplyr::mutate(
      flow = (Volume/TAUsea) - flow #from x_RegSea2Cont
    ) |>
    dplyr::select(from.Scale, to.Scale, SubCompart, Species, flow)
  
 
  return(data.frame(out))  
  
           # switch (SubCompartName,
           #         "sea" = {
           #           RegSea2Cont <- all.x_RegSea2Cont$flow[all.x_RegSea2Cont$fromSubCompart == "sea" & 
           #                                                   all.x_RegSea2Cont$fromScale == "Regional"]
  #           toVolume <- all.Volume$Volume[all.Volume$SubCompart == SubCompartName &
  #                                    all.Volume$Scale == "Continental" ]
  #           toTAUsea <- all.TAUsea$TAUsea[all.TAUsea$Scale =="Continental" ]
  #           return((toVolume/toTAUsea)-RegSea2Cont)}, 
  #         NA
  # )
  # 
}
