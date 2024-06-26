#'@title FRACa
#'@name FRaction of air in any matrix
#'@description either subFRACa or, when the main matrix, remainder after substracting subFRACs + subFRACw
#'@return 
#'@export
FRACa <- function(subFRACa, subFRACw, subFRACs, Matrix) {
  if (Matrix == "air") {
    return (1 - subFRACw - subFRACs)
  } else if (Matrix == "sediment") {
    return (0)
  } else {
    return (subFRACa)
  }
}


