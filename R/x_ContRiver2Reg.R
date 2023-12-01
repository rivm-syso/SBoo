#' @title Continental River 2 Regional
#' @name x_ContRiver2Reg
#' @param x_Runoff runoff flux
#' @param dischargeFRAC fraction of ...
#' @return River Discharge for scale Continental
#' @export
x_ContRiver2Reg <- function (ScaleName, SubCompartName, 
                           all.Runoff, RainOnFreshwater, 
                           dischargeFRAC, LakeFracRiver){
  if(ScaleName == "Continental" & SubCompartName == "river"){ #the one and only
    SumRainRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == "Continental"])
    River2sea  <- RainOnFreshwater + SumRainRunoff * (1-dischargeFRAC)
    Lake2River <- LakeFracRiver * River2sea
    
    return((RainOnFreshwater + SumRainRunoff + Lake2River) * dischargeFRAC)
    
  } else NA
}