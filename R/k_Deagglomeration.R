#' @title Fragmentation and degradation of plastics
#' @name k_Deagglomeration
#' @description Deagglomeration is the opposite of heteroagglomeration. For instance Polymer particles with an inorganic part can fragment. It is used in Plastics World. For instance for fragmentation of Tyre and Road wear particles into Tyre wear particles alone.
#' @param kdeag Falling appart or fragmenting of heteroaglomerates [s-1]
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_Deagglomeration, the combined degradation and fragmentation rate of microplastics [s-1]
#' @export


k_Deagglomeration <- function (kdeag, SubCompartName, ScaleName, parent){
  
  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Deagglomeration")))
  
  out <- parent$FromDataAndTo("k_Deagglomeration") |>
    # From
    dplyr::left_join(kdeag, by=c("SubCompart", "from.Species"="Species")) |>
    dplyr::filter(!(
      Scale %in% c("Tropic", "Moderate", "Arctic") & SubCompart %in% c("freshwatersediment", "lakesediment", "lake", "river", "agriculturalsoil", "othersoil")
    )) |>
    dplyr::mutate(
      k_Deagglomeration = kdeag
    )|>
    dplyr::filter(!is.na(k_Deagglomeration)) |>
    dplyr::select(names, k_Deagglomeration) |>
    dplyr::arrange(dplyr::across(names))
  
  return(data.frame(out))

  # if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & (SubCompartName == "freshwatersediment" | 
  #                                                            SubCompartName == "lakesediment" |
  #                                                            SubCompartName == "lake" |
  #                                                            SubCompartName == "river" |
  #                                                            SubCompartName == "agriculturalsoil"|
  #                                                            SubCompartName == "othersoil")) ){
  #   return(NA)
  # } else if (!is.null(ScaleName) && ScaleName %in% c("Regional", "Continental") && (SubCompartName == "deepocean") && (is.na(Regional_and_Continental_deepocean) || isFALSE(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE")){
  #   return(NA)
  # } else {
  #   return(kdeag)
  # }
    
}