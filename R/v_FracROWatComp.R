#' @title FracRoWatComp
#' @name  FracROWatComp
#' @description Fraction of water component compared to total water in scale for setting fraction of runoff reaching each freshwater compartment
#' @param all.landFRAC Fractions of land compartments #[-]
#' @param all.Matrix Matrix compartment the subcompartments belong to #[-]
#' @param Matrix Current matrix the function calculates for [-]
#' @param SubCompartName Subcompartment the function calculates for [-]
#' @param ScaleName The scale the function calculates for [-]
#' @return Fraction of water component of .to subcompartment #[-]
#' @export
#'


FracROWatComp <- function(landFRAC, 
                          #all.Matrix, 
                          Matrix, 
                          SubCompartName, 
                          ScaleName,
                          SpeciesName,
                          parent) {
  
  out <- ScaleName |>
    tidry::expand_grid(SubCompartName, SpeciesName) |>
    dplyr::full_join(Matrix, by="SubCompart") |>
    dplyr::full_join(landFRAC, c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    dplyr::group_by(ScaleName, Matrix) |>
    dplyr::mutate(
      FracROWatComp = dplyr::case_when(
        Matrix == 'water' & ScaleName %in% c("Regional", "Continental") ~ landFRAC / sum(landFRAC, na.rm=TRUE),
        SubCompartName == 'sea' & ScaleName %in% c("Tropic", "Moderate", "Arctic") ~ 1,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::ungroup() |>
    dplyr::filter(!is.na(FracROWatComp)) |>
    dplyr::arrange(Scale, SubCompart, Species) |>
    dplyr::select(Scale, SubCompart, FracROWatComp)
  
  return(data.frame(out)) 
  
  # if ((Matrix == "water") & (ScaleName %in% c("Regional", "Continental"))) {
  #   
  #   compFrac <- all.landFRAC$landFRAC[all.landFRAC$SubCompart == SubCompartName & all.landFRAC$Scale == ScaleName]
  #   mergeddata <- merge(all.landFRAC, all.Matrix)
  #   waterFrac <- sum(mergeddata$landFRAC[mergeddata$Matrix == "water" & mergeddata$Scale == ScaleName])
  #   return(compFrac / waterFrac)
  #   
  # } else if ((SubCompartName == "sea") & (ScaleName %in% c("Tropic", "Moderate", "Arctic"))){ 
  #   return(1)
  # }  else
  # {
  #   return(NA)
  # }
}
