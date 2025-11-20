#' @title RegSea2Cont
#' @name  x_RegSea2Cont
#' @param RiverDischarge find the flux in all. RiverDischarge fluxes
#' @param ContSea2Reg find the ContSea2Reg flux in all.(ContSea2Reg)fluxes (there is no direct relation)
#' @return River Discharge for scale Continental
#' @export
x_RegSea2Cont <- function (all.x_RiverDischarge, all.x_ContSea2Reg, Test_surface_water){
  #if isTRUE(Test_water_surface) { #connect sea to sea and deepocean to deepocean
  #  switch (SubCompartName,
   #         "sea" = {OceanCurrent},
   #         "deepocean" = {return(0)}, #fix, data in code
   #         NA
  #  )
  #}
  
  
  x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale=="Regional"]
  return(all.x_ContSea2Reg$flow + x_RiverDischarge)
}
