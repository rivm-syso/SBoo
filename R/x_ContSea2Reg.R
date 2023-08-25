#' @title ContSea2Reg
#' @name x_ContSea2Reg
#' @description there is an additional flux from continental sea to regional to compensate from the LakeFracRiver
#' @param RiverDischarge see ppt
#' @param LakeFracRiver see ppt
#' @return additional flux from continental sea to regional
#' @export
x_ContSea2Reg <- function (all.x_RiverDischarge, LakeFracRiver){
    #NB the regional river discharge determines the sea flow from continental!!
    x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale=="Regional"]
    return( ((1-LakeFracRiver)/LakeFracRiver) * x_RiverDischarge)
}
