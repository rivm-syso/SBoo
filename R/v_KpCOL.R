#'@title Kp for water colloids
#'@name KpCOL
#'@description partitioning coefficient for water colloids
#'@param  D octanol/water partitioning coefficient at neutral pH for colloids [-]
#'@param Matrix the medium, the formula is only applicable to soil and sediment
#'@export
KpCOL <- function(D, Matrix){
  
  out <- Matrix |>
    dplyr::full_join(D, by="SubCompart") |>
    dplyr::filter(Matrix == 'water') |>
    dplyr::mutate(
      KpCOL = 0.08 * D
    ) |>
    dplyr::filter(!is.na(KpCOL)) |>
    dplyr::select(SubCompart, KpCOL) |>
    dplyr::arrange(SubCompart)
  
  return(out)
  
  # if (Matrix %in% c("water")) {
  #   return(
  #     0.08*D
  #   )
  # } else return (NA)
}
