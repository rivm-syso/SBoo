#' @title Fragmentation and degradation of plastics
#' @name k_Fragmentation
#' @description Calculation of fragmentation and degradation of plastics, only used in Plastics World
#' @param kmpdeg degradation rate of plastic in certain subcompartment [s-1]
#' @param kfrag fragmentation rate of plastic in certain subcompartment
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_Fragmentation, the combined degradation and fragmentation rate of microplastics [s-1]
#' @export

k_Fragmentation <- function (kfrag, parent, SpeciesName){
  
  out <- parent$FromDataAndTo("k_Fragmentation")
  # Names of dataframe
  names <- dplyr::all_of(names(out))
  
  #Species for this process
  Scale <- out |> 
    dplyr::select(Scale) |>
    dplyr::distinct()
  
  exist_kfrag <- kfrag |>
    dplyr::cross_join(Scale) |>
    parent$states$clipStates()
  # Build
  out <- out |>
    dplyr::inner_join(exist_kfrag, by = c("from.SubCompart" = "SubCompart", "Scale", "Species"))
 
  # Calc
  out <- out |>
    dplyr::filter(!(
      Scale %in% c("Tropic", "Moderate", "Arctic") & 
        from.SubCompart %in% c("freshwatersediment", "lakesediment", "lake", "river", "agriculturalsoil", "othersoil") 
    )) |>
    dplyr::mutate(k_Fragmentation = kfrag) |>
    
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