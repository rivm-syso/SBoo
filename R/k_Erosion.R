#'@title Erosion of particles
#'@name k_Erosion
#'@description Erosion process for attached ENP species corrected for the compartment their flowing to.
#'@param relevant_depth_s relevant depth that is susceptible to erosion [m], see f_CORRsoil
#'@param penentration_depth_s assumed penetration depth [m], from Hollander et al. (2007), see f_CORRsoil
#'@param VertDistance Vertical distance of compartment [m]
#'@param EROSIONsoil Soil erosion #[m/s]
#'@param ScaleName To adjust for the absence of rivers in global scales
#'@param to.SubCompartName name of the subcompartment of the destination box of this process
#'@param Matrix compartment type considered
#'@param landFRAC fraction of land of a certain compartment 
#'@param FracROWatComp a fraction of the water that erosion/runoff can go to [-], see v_FracROWatComp.
#'@return k_Erosion Erosion of attached ENP species (P) from soil to water #[s-1]
#'@export
k_Erosion <- function(relevant_depth_s,penetration_depth_s, EROSIONsoil, VertDistance,
                      FracROWatComp,
                      ScaleName, SubCompartName,  Matrix, landFRAC, SpeciesName, parent){
  
  
  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Erosion")))

  # Build grid
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Matrix", "EROSIONsoil", "VertDistance", "landFRAC")
  dfs <- mget(args, envir = environment())
  grid_from <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  # To
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Matrix", "FracROWatComp")
  dfs <- mget(args, envir = environment())
  grid_to <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()


  # Output
  out <- parent$FromDataAndTo("k_Erosion") |>
    # From
    dplyr::left_join(grid_from, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species")) |>

    #To
    dplyr::left_join(grid_to, by=c("to.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species"), suffix=c("_from", "_to")) |>

    dplyr::mutate(
      k_Erosion = dplyr::case_when(
        Scale %in% c("Regional", "Continental") & to.SubCompart == "sea" ~ NA,
        Scale %in% c("Tropic", "Moderate", "Arctic") & to.SubCompart != "sea" ~ NA,
        TRUE ~ EROSIONsoil * f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / VertDistance * FracROWatComp  #[s-1]
      )
    ) |>
    dplyr::filter(!is.na(k_Erosion)) |>
    dplyr::select(names, k_Erosion) |>
    dplyr::arrange(dplyr::across(names))

  return(data.frame(out))
  
  
  # if (ScaleName %in% c("Regional", "Continental") & to.SubCompartName == "sea") {
  #   return(NA)
  # } 
  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "sea") {
  #   return(NA)
  # } 
  # EROSIONsoil * f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / VertDistance * to.FracROWatComp  #[s-1]
}
