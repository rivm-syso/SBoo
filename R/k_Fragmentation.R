#' @title Fragmentation and degradation of plastics
#' @name k_Fragmentation
#' @description Calculation of fragmentation and degradation of plastics, only used in Plastics World
#' @param kmpdeg degradation rate of plastic in certain subcompartment [s-1]
#' @param kfrag fragmentation rate of plastic in certain subcompartment
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_Fragmentation, the combined degradation and fragmentation rate of microplastics [s-1]
#' @export


k_Fragmentation <- function (kfrag, SubCompartName, ScaleName, 	Regional_and_Continental_deepocean, parent, SpeciesName){
  
  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Fragmentation")))
  
  # Build grid
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "kfrag")
  dfs <- mget(args, envir = environment())
  grid <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  # Output
  out <- parent$FromDataAndTo("k_Fragmentation") |>
    # From
    dplyr::left_join(grid, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species")) |>
    
    dplyr::mutate( 
      k_Fragmentation = dplyr::case_when(
        ScaleName %in% c("Tropic", "Moderate", "Arctic") & 
          SubCompartName %in% c("freshwatersediment", "lakesediment", "lake", "river", "agriculturalsoil", "othersoil") ~ NA_real_,
        ScaleName %in% c("Regional", "Continental") & SubCompartName == "deepocean" & (is.na(	Regional_and_Continental_deepocean) || isFALSE(	Regional_and_Continental_deepocean) || 	Regional_and_Continental_deepocean == "FALSE") ~ NA_real_,
        TRUE ~ kfrag
      )
    ) |>
    dplyr::filter(!is.na(k_Fragmentation)) |>
    dplyr::select(names, k_Fragmentation) |>
    dplyr::arrange(dplyr::across(names))
  
  return(data.frame(out))
  
  # if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & (SubCompartName == "freshwatersediment" | 
  #                                                            SubCompartName == "lakesediment" |
  #                                                            SubCompartName == "lake" |
  #                                                            SubCompartName == "river" |
  #                                                            SubCompartName == "agriculturalsoil"|
  #                                                            SubCompartName == "othersoil")) ){
  #   return(NA)
  # } else if (ScaleName %in% c("Regional", "Continental") && (SubCompartName == "deepocean") && (is.na(	Regional_and_Continental_deepocean) || isFALSE(	Regional_and_Continental_deepocean) || 	Regional_and_Continental_deepocean == "FALSE")){
  #   return(NA)
  # } else {
  #   return(kfrag)
  # }
    
}