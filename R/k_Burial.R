#Script for determining the burial by sediment in water

##### section to be moved to SBoo after testing 
#'@title Burial by sediment
#'@name k_Burial
#'@param VertDistance mixed depth water sediment compartment [m]
#'@param NETsedrate Net sediment accumulation rate (from the surface water above) [m/s]
#'@return k_Burial Burial from sediment [s-1]
#'@export

k_Burial <- function(VertDistance, all.NETsedrate, ScaleName, SubCompartName){
  # NETsedrate assumed identical to NETsedrate of the water column above the sediment
  waterabove <- switch (SubCompartName,
      "lakesediment" = "lake",
      "marinesediment" = {switch(ScaleName, 
                                 "Tropic" = "deepocean",
                                 "Moderate" = "deepocean",
                                 "Arctic" = "deepocean",
                                 "sea")}
        ,
      "freshwatersediment" = "river",
      NA
  )
  if (is.na(waterabove)) return (NA)
  rightRow <- which(all.NETsedrate$SubCompart == waterabove & all.NETsedrate$Scale == ScaleName)
  if (length(rightRow) != 1) return(NA)
  waterNETsedrate <- all.NETsedrate$NETsedrate[rightRow]
  waterNETsedrate / VertDistance
}
