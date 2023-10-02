#'@title FRACw
#'@name FRaction of water in any matrix
#'@description either subFRACw or, when the main matrix, remainder after substracting subFRACs + subFRACa
#'@return 
#'@export
FRACw <- function(subFRACa, subFRACw, subFRACs, Matrix){
  if (Matrix == "water") {
    return (1 - subFRACs - subFRACa)
  } else
    return (subFRACw)
}
