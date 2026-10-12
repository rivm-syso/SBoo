#'@title Volatilisation rate constant
#'@name k_Volatilisation
#'@description volatilisation process based on the two-film resistance model (Schwarzenbach et al., 1993) ISBN: 978-1-118-76723-8
#'@param FRACinf Fraction infiltration #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@returns Volatilisation rate constant [s-1]
#'@export

k_Volatilisation <- function(MTC_2w, MTC_2a, MTC_2s, Kacompw, FRorig, FRorig_spw, FRinw, Kscompw,
                             VertDistance, SpeciesName, Matrix, relevant_depth_s, penetration_depth_s,
                             SubCompartName, ScaleName, parent){ 

  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Volatilisation")))
  
  # Build grid
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Matrix", "MTC_2a", "Kacompw", "FRorig", "FRorig_spw", "Kscompw", "VertDistance", "FRinw")
  dfs <- mget(args, envir = environment())
  grid_from <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Matrix", "MTC_2w", "MTC_2s")
  dfs <- mget(args, envir = environment())
  grid_to <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  
  # Output
  out <- parent$FromDataAndTo("k_Volatilisation") |>
    # From
    dplyr::left_join(grid_from, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species")) |>
    
    #To
    dplyr::left_join(grid_to, by=c("to.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species"), suffix=c("_from", "_to")) |>
    
    dplyr::mutate( 
      k_Volatilisation = dplyr::case_when(
        Species == 'Unbound' & Matrix_from == 'water' ~
          ((MTC_2w*MTC_2a/(MTC_2w*(Kacompw*FRorig)+MTC_2a))*(Kacompw*FRorig)*FRinw)/VertDistance,
        
        Species == 'Unbound' & Matrix_from == 'soil' ~
          ((MTC_2s*MTC_2a)/(MTC_2s+MTC_2a/((Kacompw*FRorig_spw)/Kscompw))*f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s))/VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Volatilisation)) |>
    dplyr::select(names, k_Volatilisation) |>
    dplyr::arrange(dplyr::across(names))
  
  return(data.frame(out))
  
# 
#   if (SpeciesName %in% c("Molecular")) {
#     switch(Matrix,
#            "water" = { 
#              flux = (to.MTC_2w*from.MTC_2a/(to.MTC_2w*(Kacompw*from.FRorig)+from.MTC_2a))*(Kacompw*from.FRorig)*FRinw
#              return(flux/VertDistance)},
#            "soil" = { 
#              fcor <- f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s)
#              flux = (to.MTC_2s*from.MTC_2a)/(to.MTC_2s+from.MTC_2a/((Kacompw*from.FRorig_spw)/Kscompw))*fcor
#              return(flux/VertDistance)}
#     )
#   } else { 
#     return(NA)
#   }
}

