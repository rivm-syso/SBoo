#'@title FracRoWatComp
#'@name  FracROWatComp
#'@description fraction of water component compared to total water in scale for correction of runoff
#'@param all.landFRAC Fractions of land compartments #[-]
#'@param all.Matrix Matrix compartment the subcompartments belong to #[-]
#'@param Matrix Current matrix the function calculates for [-]
#'@param SubCompartName Subcompartment the function calculates for [-] 
#'@param ScaleName The scale the function calculates for [-]
#'@return Fraction of water component of .to subcompartment #[-]
#'@export
#'
#

FracROWatComp <- function(all.landFRAC, all.Matrix, Matrix, SubCompartName, ScaleName) {
  compFrac <- all.landFRAC$landFRAC[all.landFRAC$SubCompart == SubCompartName & all.landFRAC$Scale ==  ScaleName]
  all.landFrac <- as.tibble(all.landFRAC)
  all.Matrix <- as.tibble(all.Matrix)
  mergeddata <- left_join(
    x = all.landFRAC,
    y = all.Matrix,
    by = join_by(SubCompart))

  if ((Matrix == "water") & (ScaleName %in% c("Regional", "Continental"))) {
    # total landfrac of (fresh) water compartments
    waterFrac <- mergeddata |>
      filter(Matrix == "water" & Scale == ScaleName) |>
      summarise(waterFrac = sum(landFRAC, na.rm = TRUE)) |>
      pull(waterFrac)
    return(compFrac / waterFrac)
  } else {
    return(1)
  }
}