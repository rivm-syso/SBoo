#'@title Adsorption of chemical to across an interface
#'@name k_Adsorption
#'@description The adsorption rate constant (k_Adsorption) is calculated based on the ... See 2015-0161 report.
#'@param FRingas see FRingas()
#'@param FRinw see FRinw()
#'@param MTC_2sd SBvariable, see MTC_2sd()
#'@param MTC_2w SBvariable
#'@param MTC_2a SBvariable
#'@param MTC_2s SBvariable
#'@param FRorig SBvariable
#'@param Kacompw SBvariable
#'@param VertDistance data
#'@param Matrix # can use paraminherit to use description from another function! No need to copy paste same descriptions!
#'@returns The adsorption rate constant relevant for the receiving compartments soil, water or sediment [s-1]
#'@export
k_Adsorption <- function (from.FRingas, from.FRinw, from.MTC_2sd, FRorig_spw,
                          all.MTC_2w, from.MTC_2a, from.MTC_2s, FRorig, Kacompw, 
                          Kscompw, to.Matrix, VertDistance, 
                          AreaLand, AreaSea, Area,
                          from.SubCompartName, to.SubCompartName, ScaleName) {
  
  to.MTC_2w <- all.MTC_2w$MTC_2w[all.MTC_2w$SubCompart == to.SubCompartName]
  from.MTC_2w <- all.MTC_2w$MTC_2w[all.MTC_2w$SubCompart == from.SubCompartName]
  
  switch(to.Matrix, #check if to/from works, current implementation is to.Matrix. so air to water and air to soil, and water to sediment.
         
         "water" = { # air to water
           if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & (to.SubCompartName %in% c("lake", "river"))) {
             return(NA)
           } 
           GASABS = from.FRingas*(from.MTC_2w*from.MTC_2a/(from.MTC_2w*(Kacompw*FRorig)+from.MTC_2a))
           AreaFrac = Area/(AreaLand+AreaSea)
           return(GASABS/VertDistance*AreaFrac) },
         "soil" = { # air to soil
           if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & 
               (to.SubCompartName %in% c("agriculturalsoil", "othersoil"))) {
             return(NA)
           } 
           GASABS = from.FRingas*(from.MTC_2s*from.MTC_2a)/(from.MTC_2s*(Kacompw*FRorig_spw)/Kscompw+from.MTC_2a)
           AreaFrac = Area/(AreaLand+AreaSea)
           return(GASABS/VertDistance*AreaFrac) },
         "sediment" = { # water to sediment
           if (ScaleName %in% c("Arctic", "Moderate", "Tropic") & from.SubCompartName != "sea") {
             return(NA)
           }
           ADSORB = (from.MTC_2sd*to.MTC_2w)/(from.MTC_2sd+to.MTC_2w)*from.FRinw
           return(ADSORB/VertDistance) }, 
         return(NA)
  )
  
}
