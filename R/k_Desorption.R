#'@title Desorbtion of particles 
#'@name Desorption
#'@param Ksdcompw
#'@param MTC_2w i.e. sediment
#'@param MTC_2sd i.e. water
#'@param VertDistance
#'@return Desorption rate constant from sediment to water [s-1]
#'@export
k_Desorption <- function (Ksdcompw, MTC_2w, MTC_2sd, VertDistance,
                          SpeciesName, to.SubCompartName, ScaleName) {
  
  switch (SpeciesName,
    "molcular" = {
      if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName == "sea") {
        return(NA)
      }
      ( (MTC_2sd*MTC_2w)/(MTC_2sd + MTC_2w)/Ksdcompw ) /
        VertDistance
    },
    return(NA)
  )
  

  
}

