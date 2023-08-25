#'@title PARTIAL MASS TRANSFER COEFFICIENT water to sediment
#'@name MTC_2sd
#'@description 
#'@param kwsd.water Constant water side pMTC to sediment [m/s-1]
#'@param from.Matrix
#'@return 
#'@export
MTC_2sd <- function(kwsd.water, Matrix){
  ret <- switch(Matrix,
    "water" =  kwsd.water,
    NA
  )
}
