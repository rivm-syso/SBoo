#' @title RegSea2Cont
#' @name  x_RegSea2Cont
#' @param RiverDischarge find the flux in all. RiverDischarge fluxes
#' @param ContSea2Reg find the ContSea2Reg flux in all.(ContSea2Reg)fluxes (there is no direct relation)
#' @return River Discharge for scale Continental
#' @export
x_RegSea2Cont <- function (x_RiverDischarge, x_ContSea2Reg, parent){
  
  out <- parent$FromDataAndTo("x_RegSea2Cont") |>
    dplyr::left_join(x_RiverDischarge, by=c("from.Scale" = "Scale", "SubCompart" = "to.SubCompart")) |>
    dplyr::left_join(x_ContSea2Reg, by=c("to.Scale" = "from.Scale", "SubCompart" = "SubCompart")) |>
    dplyr::mutate(
      x_RegSea2Cont = x_ContSea2Reg + x_RiverDischarge 
    ) |>
    dplyr::select(from.Scale, to.Scale, SubCompart, Species, x_RegSea2Cont) |>
    dplyr::filter(!is.na(x_RegSea2Cont))
    
  return(data.frame(out))
  
  # x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale=="Regional"]
  # return(all.x_ContSea2Reg$flow + x_RiverDischarge)
}
