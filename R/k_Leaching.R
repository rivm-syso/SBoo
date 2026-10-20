#'@title Leaching of particles due to infiltration into deeper soil layers
#'@name k_Leaching
#'@param FRACinf Fraction infiltration #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance depth compartment #[m]
#'@param SpeciesName considered species 
#'@param penetration_depth_s penetration depth for correction factor, based on Hollander et al. (2007) https://doi.org/10.1016/j.envpol.2006.09.018 [m]
#'@param Kscompw dimensionless partitioning coefficient soil-water [-]
#'@description Note: The fraction infiltration (FRACinf) is rescaled to the rest of inflow from rain based on the fraction rain going to runoff (FRACrun). 
#'@return Leaching rate constant of compound from  soil out of the system (e.g. no deep(er) soil layer considered). [s-1]
#'@export
k_Leaching <- function(FRACinf,FRACrun, RAINrate, VertDistance, SpeciesName, 
                       penetration_depth_s, Kscompw, parent){ #k_ Leaching
  
  # Attached does not leach??
  
  out <- parent$FromDataAndTo("k_Leaching")
  # Names of dataframe
  names <- dplyr::all_of(names(out))
  
  # Build 
  args <- c("out", "SpeciesName", "FRACinf", "RAINrate", "VertDistance", "Kscompw")
  dfs <- mget(args, envir = environment())
  out <- purrr::reduce(dfs, ~ dplyr::left_join(.x, .y))
  
  # Output
  out <- out |>
    
    dplyr::mutate( 
      
      penetration_depth_s = dplyr::if_else(SpeciesName == "Molecular", 0.1, penetration_depth_s),
      CORRleach = f_CORRsoil(VertDistance, relevant_depth_s=0.5, penetration_depth_s),
      
      #If due to extreme inputs, FRACinf and FRACrun are both too large, scale FRACinf to mass conserve rain
      FRACinf = dplyr::case_when(
        (FRACinf + FRACrun) > 1 ~ 1-FRACrun,
        TRUE ~ FRACinf
      ),
      k_Leaching = dplyr::case_when(
        SpeciesName %in%  c("Aggregated", "Nanoparticle") ~ (FRACinf * RAINrate * CORRleach) / VertDistance,
        SpeciesName %in% c("Molecular") ~ FRACinf*RAINrate/Kscompw*CORRleach/VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Leaching)) |>
    dplyr::select(all_of(names), k_Leaching) |>
    dplyr::arrange(dplyr::across(names))
  
  return(data.frame(out))
  
  # if(SpeciesName %in%  c("Aggregated", "Nanoparticle") ){
  #   #Correction factor depth dependent soil concentration
  #   CORRleach <- f_CORRsoil(VertDistance, relevant_depth_s=0.5, penetration_depth_s) 
  #   #If due to extreme inputs, FRACinf and FRACrun are both too large, scale FRACinf to mass conserve rain
  #   if((FRACinf + FRACrun) > 1) {FRACinf = 1-FRACrun}
  #   #Leaching of aggregated (A) and free (S) enp species from natural soil [s-1]
  #   return( (FRACinf * RAINrate * CORRleach) / VertDistance )
  # } else if( SpeciesName %in% c("Molecular")) {
  #   CORRleach <- f_CORRsoil(VertDistance, relevant_depth_s=0.5, penetration_depth_s=0.1)
  #   #If due to extreme inputs, FRACinf and FRACrun are both too large, scale FRACinf to conserve rain
  #   if((FRACinf + FRACrun) > 1) {FRACinf = 1-FRACrun}
  #   return(FRACinf*RAINrate/Kscompw*CORRleach/VertDistance)
  # } else return(NA)
  
}
