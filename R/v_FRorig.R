#'@title Fraction of orginigal species - air, soil, sediment and water
#'@name FRorig
#'@description Fraction original species in air, soil, sediment or water
#'@param pH pH of soil, sediment, water or aerosol water
#'@param pKa Dissociation constant of (conjugated) acid (default = 7)
#'@param ChemClass Class of chemical, in this case Acid, Base or other
#'@param Matrix type of compartment considered
#'@return FRorig
#'@export
FRorig <- function(ChemClass, Matrix,pH, pKa){
  
  out <- Matrix |> 
    tidyr::expand_grid(pKa) |>
    dplyr::full_join(pH, by="SubCompart") |>
    dplyr::mutate(
      FRorig = dplyr::case_when(
        (Matrix == 'soil' | Matrix == 'sediment') & ChemClass == 'acid' & !is.na(pKa) ~ 1/(1+10^(pH-0.6-pKa)),
        (Matrix == 'soil' | Matrix == 'sediment') & ChemClass == 'base' & !is.na(pKa) ~ 1/(1+10^(pKa-4.5)),
        (Matrix == 'water' | Matrix == 'air') & ChemClass == 'acid' & !is.na(pKa) ~ 1/(1+10^(pH-pKa)),
        (Matrix == 'water' | Matrix == 'air') & ChemClass == 'base' & !is.na(pKa) ~ 1/(1+10^(pKa-pH)),
        (Matrix == 'water' | Matrix == 'sediment' | Matrix == 'air' | Matrix == 'soil') ~ 1,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(FRorig)) |>
    dplyr::select(SubCompart, FRorig)
  
  return(data.frame(out))       
}
