#'@title Desorbtion of molecular species from sediment to water
#'@name k_Desorption
#'@param ...
#'@return Desorption rate constant from sediment to water [s-1]
#'@export
k_Desorption <- function (Ksdcompw, MTC_2w, to.MTC_2sd, VertDistance,
                          SpeciesName, to.SubCompartName, ScaleName) {
  
  switch (SpeciesName,
    "Molecular" = {
      # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName == "sea") {
      #   return(NA)
      # }
      ( (to.MTC_2sd*MTC_2w)/(to.MTC_2sd + MTC_2w)/Ksdcompw ) /
        VertDistance
    },
    return(NA)
  )
  
}

