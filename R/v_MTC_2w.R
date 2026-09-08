#'@title PARTIAL MASS TRANSFER COEFFICIENT air and sediment to water
#'@name MTC_2w
#'@description mass transfer coefficient to water
#'@param WINDspeed Windspeed in compartment/scale [m.s-1]
#'@param MW Molecular weight of compound kg [g.mol-1]
#'@param kwsd.sed Constant sediment transfer rate to water [m/s-1]
#'@param from.Matrix Matrix/compartment from which the relevant process is taking place
#'@return MTC_2w
#'@export
MTC_2w <- function(WINDspeed, MW, kwsd.sed, Matrix, parent, ScaleName, SpeciesName){
  
  out <- ScaleName |>
    tidyr::expand_grid(Matrix, SpeciesName) |>
    parent$states$clipStates() |>
    dplyr::filter(Matrix %in% c("air", "sediment")) |>
    dplyr::full_join(WINDspeed) |>
    dplyr::mutate(
      MTC_2w = dplyr::case_when(
        Matrix == 'air' ~ 0.01*(0.3+0.2*WINDspeed)*((0.018/MW)^(0.67*0.5)),
        Matrix == 'sediment' ~  kwsd.sed,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::arrange(Scale, SubCompart) |>
    dplyr::select(Scale, SubCompart, MTC_2w)
  
  return(data.frame(out))
}
