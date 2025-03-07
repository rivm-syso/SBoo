#' @title Continental surface water to Regional [s-1]
#' @name x_ContRiver2Reg
#' @param ScaleName Name of the relevant scale
#' @param SubCompartName Name of the relevant sub-compartment
#' @param Runoff Runoff flow from soil to river [m3.s-1]
#' @param RainOnFreshwater Water flow of rain directly on lake/river [m3.s-1]
#' @param dischargeFRAC Fraction discharge regional fresh water to continental scale [-]
#' @param LakeFracRiver Fraction of river discharge that originates from lake [-]
#' @return River Discharge for scale Continental [s-1]
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