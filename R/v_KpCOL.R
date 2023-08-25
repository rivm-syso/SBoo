#'@title Kp for water colloids
#'@name KpCOL
#'@description subcompartment/water PARTITION COEFFICIENT
#'@param D
#'@param Matrix the medium, the formula is only applicable to soil and sediment
#'@export
KpCOL <- function(D, Matrix){
  if (Matrix %in% c("water")) {
    return(
      0.08*D
    )
  } else return (NA)
}
