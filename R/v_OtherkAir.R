#' @title v_OtherkAir
#' @name v_OtherkAir
#' @param all.kaas sneaky way to obtain other loss processes from air
#' @export
OtherkAir <- function (all.kaas){
  OnlyMolecular <- all.kaas$fromSpecies == "Unbound" # Unbound and not Molecular as this is not SpeciesName
  OnlyCalculated <- all.kaas$process != "LoadKaas" # Not sure if this is needed here, but at moment all.kaas consists also of LoadKaas which are doubled with the calculated K's
  
  # TODO incorporate OnlyCalculated in selection below
  
  return(data.frame(
      Scale = all.kaas$fromScale[all.kaas$fromSubCompart == "air" & OnlyMolecular] ,
      Species = all.kaas$fromSpecies[all.kaas$fromSubCompart == "air" & OnlyMolecular],
      OtherkAir = all.kaas$k[all.kaas$fromSubCompart == "air" & OnlyMolecular]
    ))
  
  
  
}
