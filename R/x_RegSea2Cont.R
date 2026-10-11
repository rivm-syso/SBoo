#' @title RegSea2Cont
#' @name  x_RegSea2Cont
#' @param RiverDischarge find the flux in all. RiverDischarge fluxes
#' @param ContSea2Reg find the ContSea2Reg flux in all.(ContSea2Reg)fluxes (there is no direct relation)
#' @return River Discharge for scale Continental
#' @export
x_RegSea2Cont <- function (x_RiverDischarge, x_ContSea2Reg, parent){
  
  out <- parent$FromDataAndTo("x_RegSea2Cont") |> #
  #   from.Scale    to.Scale SubCompart Species
  # 2   Regional Continental  deepocean Unbound to deepocean, really? for now filter after next
    dplyr::left_join(x_RiverDischarge, by=c("from.Scale" = "Scale", "SubCompart" = "to.SubCompart", "Species")) |>
    dplyr::select(-from.SubCompart) |>
    dplyr::filter(!is.na(flow)) |>
    dplyr::left_join(x_ContSea2Reg, by=c("to.Scale" = "from.Scale", "SubCompart", "Species")) |>
    dplyr::select(-to.Scale.y) |>
    dplyr::mutate(
      x_RegSea2Cont = x_ContSea2Reg + flow # from x_RiverDischarge 
    ) |>
    dplyr::select(from.Scale, to.Scale, SubCompart, Species, x_RegSea2Cont) |>
    dplyr::rename(flow = x_RegSea2Cont) 
    
  return(data.frame(out))
  
  # x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale=="Regional"]
  # return(all.x_ContSea2Reg$flow + x_RiverDischarge)
}
