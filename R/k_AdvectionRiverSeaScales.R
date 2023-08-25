#' @title k_AdvectionRiverSeaScales
#' @name k_AdvectionRiverSeaScales
#' @param from.SubCompartName dimension dependency in function
#' @param from.ScaleName dimension dependency in function
#' @param to.ScaleName dimension dependency in function
#' @param from.Volume Volume [m3] of the from box
#' @param to.Volume Volume [m3] of the from box
#' @param from.TAUsea refresh rate
#' @param to.TAUsea dito for the "to" box
#' @param ContRiver2Reg riverflow from continental to regional
#' @param ContSea2Reg seaflow from continental to regional
#' @param RegSea2Cont seaflow from regional to continental
#' @param OceanCurrent global constant
#' @param ContinentalInModerate logical; continental can be in Moderate or in Tropic
#' @return Advection seas and deepoceans, excluding vertical mixing sea <-> deepocean
#' @export
k_AdvectionRiverSeaScales <- function (
                     from.SubCompartName,
                     from.ScaleName,
                     to.ScaleName,
                     from.Volume, to.Volume,
                     ContRiver2Reg,
                     ContSea2Reg,
                     RegSea2Cont,
                     OceanCurrent,
                     from.TAUsea, to.TAUsea,
                     ContinentalInModerate
  ) {
    #Dims <- list(...) No longer needed; use dimension Name

      #1) between Regional and Continental scale 
        #1a river 2 river
        #1b sea 2 sea
    if (from.SubCompartName == "river" & from.ScaleName == "Continental") {
      return(ContRiver2Reg / from.Volume)
    }
    
    if (from.ScaleName == "Regional") {
      if (to.ScaleName == "Continental") return(RegSea2Cont / from.Volume)

    } else if (from.ScaleName == "Continental") {
      if (to.ScaleName == "Regional") {
        return(ContSea2Reg / from.Volume)
      } else {
        if ( (to.ScaleName == "Moderate" & ContinentalInModerate) |
             (to.ScaleName == "Tropic" & !ContinentalInModerate) ) {
          
          return((from.Volume / from.TAUsea - RegSea2Cont) / from.Volume) #yes, flow from regional / continental volume
        } 
      }
    } else { #global scales
      
      #2) Sea continental within it's scale
      if (from.SubCompartName == "sea"){
        if ( (to.ScaleName == "Continental" & from.ScaleName == "Moderate" & ContinentalInModerate) |
             (to.ScaleName == "Continental" & from.ScaleName == "Tropic" & !ContinentalInModerate) ) {
          # same flow as the opposite flow
          FlowFromContinent <-  to.Volume / to.TAUsea - RegSea2Cont 
          return(FlowFromContinent / from.Volume)
        } else {
          
      #3) Seas between global scales
          # OceanCurrent is circular flow (NB only sea-sea here!): 
          #   moderate deepocean <- moderate sea <- tropic sea <- tropic deepocean <- moderate deepocean and
          #   moderate deepocean -> moderate sea -> arctic sea -> arctic deepocean -> moderate deepocean 
          # this function implement the sea-sea part only; more is in k_AdvectionSeaOcean
          # the opposite way for arctic and tropic is counterflow, so k <- 0 (or absent)
          if ( (to.ScaleName == "Moderate" & from.ScaleName == "Tropic") |
               (to.ScaleName == "Arctic" & from.ScaleName == "Moderate") ) {
            return(OceanCurrent / from.Volume)
          } else {
            return(NA)
          }
        }
        
      } else {
        
      #4 Oceans between global scales; opposite to sea current
        if ( (from.ScaleName == "Moderate" & to.ScaleName == "Tropic") |
             (from.ScaleName == "Arctic" & to.ScaleName == "Moderate") ) {
          return(OceanCurrent / from.Volume)
        } else {
          return(NA)
        }
        
      }
    }
    
  #all other cases 
  return(NA)
}