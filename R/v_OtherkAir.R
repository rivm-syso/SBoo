#' @title v_OtherkAir
#' @name v_OtherkAir
#' @param all.kaas sneaky way to obtain other loss processes from air
#' @export
OtherkAir <- function (all.kaas){
  OnlyMolecular <- all.kaas$fromSpecies == "Molecular"
  return(data.frame(
      Scale = all.kaas$fromScale[all.kaas$fromSubCompart == "air" & OnlyMolecular] ,
      Species = all.kaas$fromSpecies[all.kaas$fromSubCompart == "air" & OnlyMolecular],
      OtherkAir = all.kaas$k[all.kaas$fromSubCompart == "air" & OnlyMolecular]
    ))
}
