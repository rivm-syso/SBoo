#' @title Volume of the SubCompartment
#' @name Volume
#' @param VertDistance hight or depth for subcompartments below the horizon
#' @param Area in m2
#' @param FRACcldw fraction of cloudwater
#' @return Volume
#' @export
Volume <- function (VertDistance, Area, FRACcldw, SubCompartName){
  
  if(SubCompartName == "air"){
    VertDistance * Area * (1-FRACcldw)
  } else if(SubCompartName == "cloudwater") {
    VertDistance * Area * FRACcldw
  } else 
    VertDistance * Area 
}
