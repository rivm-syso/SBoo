#'@title PARTIAL MASS TRANSFER COEFFICIENT water to sediment
#'@name MTC_2sd
#'@description partial mass transfer coefficient water to sediment
#'@param kwsd.water Constant water transfer to sediment [m s-1]
#'@param Matrix type of compartment
#'@param SubCompartName name of subcompartment
#'@return MTC_2sd
#'@export
MTC_2sd <- function(kwsd.water, Matrix, SubCompartName){
  
  out <- dplyr::full_join(SubCompartName, Matrix) |>
    dplyr::filter(SubCompart != 'cloudwater' & Matrix == 'water') |>
    mutate(
      MTC_2sd = kwsd.water
    ) |>
    filter(!is.na(MTC_2sd)) |>
    arrange(SubCompart) |>
    select(SubCompart, MTC_2sd)
  
  return(out)
}
