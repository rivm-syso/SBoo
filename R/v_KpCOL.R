#'@title Kp for water colloids
#'@name KpCOL
#'@description partitioning coefficient for water colloids
#'@param  D octanol/water partitioning coefficient at neutral pH for colloids [-]
#'@param Matrix the medium, the formula is only applicable to soil and sediment
#'@export
KpCOL <- function(D, Matrix){
  
  out <- Matrix |>
    full_join(D, by="SubCompart") |>
    filter(Matrix == 'water') |>
    mutate(
      KpCOL = 0.08 * D
    ) |>
    filter(!is.na(KpCOL)) |>
    select(SubCompart, KpCOL) |>
    arrange(SubCompart)
  
  return(out)
  
  # if (Matrix %in% c("water")) {
  #   return(
  #     0.08*D
  #   )
  # } else return (NA)
}
