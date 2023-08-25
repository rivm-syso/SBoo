#'@title PARTIAL MASS TRANSFER COEFFICIENT air to soil
#'@name MTC_2s
#'@description 
#'@param Mackay1 constant described by Mackay (2001) [m/s] # data in m/d (units conversion)
#'@param Mackay2 constant described by Mackay (2001) [-]
#'@param from.Matrix
#'@return 
#'@export
MTC_2s <- function(Mackay1, Mackay2, Matrix){
  ret <- switch(Matrix,
    "air" =  Mackay1/Mackay2,
    NA
  )
}
