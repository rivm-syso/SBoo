#' @title k_AdvectionWaters
#' @name k_AdvectionWaters
#' @param WATERflow.w0R.w2R 
#' @param WATERflow.w1R.w0R 
#' @param Volume of box, SB variable
#' @param LakeOutflow flow
#' @param RiverDischarge SB variable
#' @return Advection lake - river - to sea; only local scales
#' @export
k_AdvectionWaters <- function (WATERflow.w0R.w2R,
                               WATERflow.w1R.w0R, 
                     Volume,
                     LakeOutflow,
                     RiverDischarge) {

    stopifnot(Volume > 0)
    
    if (Dims$fromSubCompart == "lake") {
      if (Dims$toSubCompart == "sea") {
        return(WATERflow.w0R.w2R / Volume)
      } else if (Dims$toSubCompart == "river") {
        return(LakeOutflow / Volume)
      } 
    }

    if (Dims$fromSubCompart == "river") {
      if (Dims$toSubCompart == "sea") {
        return(RiverDischarge / Volume)
      } else if (Dims$toSubCompart == "lake") {
        return(WATERflow.w1R.w0R / Volume)
      }
    }
    
    #all other cases 
    stop(paste ("unknown k_AdvectionWaters from", Dims$SubCompart,
            Dims$Scale, "to",  Dims$toSubCompart, Dims$toScale))

}