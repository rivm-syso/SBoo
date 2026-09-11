#' @title Ocean currents flowing to Moderate scale (w2 and w3) 
#' @description 
#' Advection of continental sea, Tropic and Arctic sea compartments to sea and deepocean at Moderate scale.
#' @param Volume Volume of compartment [m3]
#' @param TAUsea Residence time of water in sea - scale variable [s]
#' @param all.x_RegSea2Cont All flows from regional to continental seawater [s-1]
#' @param OceanCurrent Global ocean circulation current [m3.s-1] 
#' @param SubCompartName Name of the subcompartment of the box at hand
#' @param ScaleName Name of the scale of the box at hand
#' @param Remove_global If this variable is TRUE, the global scales (Arctic, Moderate and Tropic) are removed
#' @return Water flow to Moderate scale surface and deepocean waters [m3 s-1]
#' @export
#' 
x_ToModerateWater <- function (Volume, TAUsea, parent, x_RegSea2Cont, OceanCurrent, Remove_global) {
  
  x_RegSea2Cont <- dplyr::filter(x_RegSea2Cont, SubCompart == 'sea', from.Scale == 'Regional') |> dplyr::pull(x_RegSea2Cont)
  
  out <- parent$FromDataAndTo("x_ToModerateWater") |>
    dplyr::left_join(Volume, by=c("from.Scale" = "Scale", "SubCompart" = "SubCompart")) |> 
    dplyr::left_join(TAUsea, by=c("from.Scale" = "Scale")) |>
    dplyr::mutate(
      x_ToModerateWater = dplyr::case_when(
        SubCompart == 'sea' & from.Scale == "Arctic" ~ 0,
        
        SubCompart == 'deepocean' & from.Scale == "Arctic" ~ OceanCurrent,
        
        SubCompart == 'sea' & from.Scale == "Tropic" ~ OceanCurrent,
        
        SubCompart == 'deepocean' & from.Scale == "Tropic" ~ 0,
        
        SubCompart == 'sea' & from.Scale == 'Continental' & 
          ((!is.na(Remove_global) && (!isTRUE(Remove_global) || Remove_global == "FALSE")))  ~ 
          (Volume/TAUsea) - x_RegSea2Cont,
        
        from.Scale == 'Continental' & SubCompart == 'sea' & 
          ((!is.na(Remove_global) && (isTRUE(Remove_global) || Remove_global == "TRUE"))) ~
          NA,
        TRUE ~ NA
      )
    ) |>
    dplyr::select(from.Scale, to.Scale, SubCompart, Species, x_ToModerateWater)
  
  return(data.frame(out))
  
  # switch(ScaleName,
  #        "Tropic" = {
  #          switch (SubCompartName,
  #                  "sea" = {OceanCurrent},
  #                  "deepocean" = {return(0)}, #fix, data in code
  #                  NA
  #          )},
  #        "Arctic" = {
  #          switch (SubCompartName,
  #                  "sea" = {return(0)}, #fix, data in code
  #                  "deepocean" = {OceanCurrent},
  #                  NA
  #          )},
  #        "Continental" = {
  #          switch (SubCompartName,
  #                  "sea" = {
  #                    if ((!is.na(Remove_global) && (isTRUE(Remove_global) || Remove_global == "TRUE"))) {
  #                      return(NA) } 
  #                    RegSea2Cont <- all.x_RegSea2Cont$flow[all.x_RegSea2Cont$fromSubCompart == "sea" & 
  #                                                            all.x_RegSea2Cont$fromScale == "Regional"]
  #                    return((Volume/TAUsea)-RegSea2Cont)}, 
  #                  NA
  #          )},
  #        return(NA)
  # )
  
}
