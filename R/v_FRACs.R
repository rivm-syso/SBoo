#'@title FRACs
#'@name FRaction of solid in any matrix
#'@description either subFRACs or, when the main matrix, remainder after substracting subFRACa + subFRACw
#'@return 
#'@export
FRACs <- function(FRACa, FRACw, FRACs, Matrix){
  if (Matrix %in% c("soil", "sediment")) {
    if (Matrix == "sediment") FRACa <- 0 #not in data
    return (1 - FRACw - FRACa)
  } else
    return (FRACs)
}
