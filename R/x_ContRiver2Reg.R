#' @title Continental surface water to Regional
#' @name x_ContRiver2Reg
#' @param Runoff runoff flow from soil to river
#' @param dischargeFRAC fraction of discharge that is going from continental to regional surface water
#' @param RainOnFreshwater 
#' @param LakeFracRiver description
#' @param ScaleName Name of the relevant scale
#' @param SubCompartName Name of the relevant sub-compartment
#' @return River Discharge for scale Continental
#' @export
x_ContRiver2Reg <- function (ScaleName, SubCompartName, 
                             all.Runoff, RainOnFreshwater, 
                             dischargeFRAC, LakeFracRiver){
  switch (ScaleName,
          "Continental" = {
            switch (SubCompartName,
                    "river" = {            SumRainRunoff <- sum(all.Runoff$Runoff[all.Runoff$Scale == "Continental"])
                    River2sea  <- RainOnFreshwater + SumRainRunoff * (1-dischargeFRAC)
                    Lake2River <- LakeFracRiver * River2sea
                    
                    return((RainOnFreshwater + SumRainRunoff + Lake2River) * dischargeFRAC)
                    },
                    return(NA)
            )
          },
          return(NA)
            )

}