#' @title ContSea2Moder
#' @name  x_ContSea2Moder
#' @param Volume in m3
#' @param TAUsea average residence time 
#' @param x_RegSea2Cont flux from regional sea to continental sea
#' @return water flux; circulation (V/Tau) - sea flux Regional 2 Continental (compensating flux originating from rainwater ??)
#' @export
x_ContSea2Moder <- function (Volume, TAUsea, all.x_RegSea2Cont, ScaleName, SubCompartName){
  RegSea2Cont <- all.x_RegSea2Cont$flow[all.x_RegSea2Cont$fromSubCompart == "sea" & 
                                          all.x_RegSea2Cont$fromScale == "Regional"]
  if(ScaleName == "Continental" & SubCompartName == "sea"){
    return((Volume/TAUsea)-RegSea2Cont)
  } else NA
}
