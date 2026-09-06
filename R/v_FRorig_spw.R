#'@title Fraction of original species - soil pore water
#'@name FRorig_spw
#'@description Fraction original species in soil pore water
#'@param pH pH of soil
#'@param pKa Dissociation constant of (conjugated) acid (default = 7)
#'@param ChemClass Class of chemical, in this case Acid, Base or other
#'@param Matrix type of compartment considered
#'@export
FRorig_spw <- function(ChemClass, Matrix, pH, pKa){
  out <- Matrix |> 
    full_join(pH, by="SubCompart") |>
    mutate(
      FRorig_spw = dplyr::case_when(
        Matrix == 'soil' & ChemClass == 'acid' & !is.na(pKa) ~ 1/(1+10^(pH-pKa)),
        Matrix == 'soil' & ChemClass == 'base' & !is.na(pKa) ~ 1/(1+10^(pKa-pH)),
        Matrix == 'soil' ~ 1,
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(FRorig_spw)) |>
    select(SubCompart, FRorig_spw)
  
  return(out)                              
}
