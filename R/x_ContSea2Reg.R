#' @title ContSea2Reg
#' @name x_ContSea2Reg
#' @description there is an additional flux from continental sea to regional to compensate from the LakeFracRiver
#' @param RiverDischarge river discharge rate, see x_RiverDischarge [s-1]
#' @return additional flux from continental sea to regional
#' @export
x_ContSea2Reg <- function (all.x_RiverDischarge){
    #NB the regional river discharge determines the sea flow from continental!!
    x_RiverDischarge <- all.x_RiverDischarge$flow[all.x_RiverDischarge$fromScale=="Regional"]
    
    return( (10-1) * x_RiverDischarge) # 10-1 based on SB4.01, maybe change to TAU based approach!
}
