#' @title v_OtherkAir
#' @name v_OtherkAir
#' @param all.kaas sneaky way to obtain other loss processes from air
#' @export
OtherkAir <- function (all.kaas, ScaleName, SpeciesName){
  #is only for Molecular
  if (SpeciesName != "Molecular") return (NA)
  #else
  #TODO clean fromScale comparison with ScaleName; which is strictly speaking not 1:1
  return(sum(all.kaas$k[all.kaas$fromScale == ScaleName & all.kaas$fromSpecies == "Unbound" & all.kaas$fromSubCompart == "air"]))

}
