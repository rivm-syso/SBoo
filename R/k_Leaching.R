#'@title Leaching of particles due to infiltration into deeper soil layers
#'@name k_Leaching
#'@param FRACinf Fraction infiltration #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@return Leaching of aggregated (A) or free (S) enp species from natural soil #[s-1]
#'@export

k_Leaching <- function(FRACinf, RAINrate, VertDistance, SpeciesName, 
                       penetration_depth_s, Kscompw){ #k_ Leaching
  
  if(SpeciesName %in%  c("Aggregated", "Nanoparticle") ){
    #Correction factor depth dependent soil concentration
    CORRleach <- f_CORRsoiltop(VertDistance, relevant_depth_s=0.5, penetration_depth_s) 
    
    #Leaching of aggregated (A) and free (S) enp species from natural soil [s-1]
    return( (FRACinf * RAINrate * CORRleach) / VertDistance )
  } else if( SpeciesName %in% c("Molecular")) {
    CORRleach <- f_CORRsoiltop(VertDistance, relevant_depth_s=0.5, penetration_depth_s) 
    return(FRACinf*RAINrate/Kscompw*CORRleach/VertDistance)
  } else return(NA)
  
}