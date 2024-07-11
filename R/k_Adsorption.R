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
#'@param to.Area Area of the receiving compartent
#'@param VertDistance data
#'@param Matrix # can use paraminherit to use description from another function! No need to copy paste same descriptions!
#'@returns The adsorption rate constant relevant for the receiving compartments soil, water or sediment [s-1]
#'@export
k_Adsorption <- function (FRingas, FRinw, from.MTC_2sd, to.FRorig_spw,
                          to.MTC_2w, from.MTC_2w, to.MTC_2a, from.MTC_2s, to.FRorig, Kacompw, 
                          to.Kscompw, to.Matrix, VertDistance, 
                          AreaLand, AreaSea, to.Area,
                          from.SubCompartName, to.SubCompartName, ScaleName, Test) {
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & from.SubCompartName == "sea") {
    return(NA)
  }
  switch(to.Matrix,
         
         "water" = { # air to water
           if (ScaleName %in% c("Regional", "Continental")){
             if (as.character(Test) == "TRUE"){
               to.FRorig <- 0.431358608116332
             } 
           }
           GASABS = FRingas*(from.MTC_2w*to.MTC_2a/(from.MTC_2w*(Kacompw*to.FRorig)+to.MTC_2a))
           AreaFrac = to.Area/(AreaLand+AreaSea)
           return(GASABS/VertDistance*AreaFrac) },
         "soil" = { # air to soil
           if (as.character(Test) == "TRUE"){
             to.FRorig_spw <- 0.00752866498600889
           } 
           GASABS = FRingas*(from.MTC_2s*to.MTC_2a)/(from.MTC_2s*(Kacompw*to.FRorig_spw)/to.Kscompw+to.MTC_2a)
           AreaFrac = to.Area/(AreaLand+AreaSea)
           return(GASABS/VertDistance*AreaFrac) },
         "sediment" = { # water to sediment
           ADSORB = (from.MTC_2sd*to.MTC_2w)/(from.MTC_2sd+to.MTC_2w)*FRinw
           return(ADSORB/VertDistance) }, 
         return(NA)
  )
  
}
