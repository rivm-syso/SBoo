#'@title FRACs
#'@name FRaction of solid in any matrix
#'@description either subFRACs or, when the main matrix, remainder after substracting subFRACa + subFRACw
#'@return 
#'@export
FRACs <- function(subFRACa, subFRACw, subFRACs, Matrix){
  if (Matrix %in% c("soil", "sediment")) {
    if (Matrix == "sediment") subFRACa <- 0
    return (1 - subFRACw - subFRACa)
  } else
    return (subFRACs)
}

