#'@title Volatilisation rate constant
#'@name k_Volatilisation
#'@param FRACinf Fraction infiltration #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@returns Volatilisation rate constant [s-1]
#'@export

k_Volatilisation <- function(to.MTC_2w, MTC_2a, to.MTC_2s, Kacompw, FRorig, FRinw, Kscompw,
                             VertDistance, SpeciesName, Matrix, relevant_depth_s, penetration_depth_s){ 
  
  if (SpeciesName %in% c("Molecular")) {
    switch(Matrix,
           "water" = { 
             flux = (to.MTC_2w*MTC_2a/(to.MTC_2w*(Kacompw*FRorig)+MTC_2a))*(Kacompw*FRorig)*FRinw
             return(flux/VertDistance)},
           "soil" = { 
             flux = (to.MTC_2s*MTC_2a)/(to.MTC_2s+MTC_2a/((Kacompw*FRorig)/Kscompw))*f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s)
             return(flux/VertDistance)}
    )
    
    
  } else { 
    return(NA)
  }
}

