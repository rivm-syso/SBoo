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


FracROWatComp <- function(all.landFRAC, 
                          #all.Matrix, 
                          Matrix, 
                          SubCompartName, 
                          ScaleName, 
                          parent) {
  
  out <- ScaleName |>
    expand_grid(SubCompartName) |>
    full_join(Matrix, by="SubCompart") |>
    full_join(all.landFRAC, c("Scale", "SubCompart")) |>
    mutate(Species ='Unbound') |>
    parent$states$clipStates() |>
    group_by(ScaleName, Matrix) |>
    mutate(
      FracRoWatComp = dplyr::case_when(
        Matrix == 'water' & ScaleName %in% c("Regional", "Continental") ~ landFRAC / sum(landFRAC, na.rm=TRUE),
        SubCompartName == 'sea' & ScaleName %in% c("Tropic", "Moderate", "Arctic") ~ 1,
        TRUE ~ NA_real_
      )
    ) |>
    ungroup() |>
    filter(!is.na(FracRoWatComp)) |>
    arrange(Scale, SubCompart, Species) |>
    select(Scale, SubCompart, Species, FracRoWatComp)
  
  return(out)  
  
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
