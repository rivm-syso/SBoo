#' @title Aerosol solids - air Partition coefficient
#' @name Kaers
#' @description Fraction of some fraction
#' @param Corg
#' @param RHOsolid
#' @param Kow
#' @return Kaers
#' @export
Kaers <- function (Kaw25,Kow, Corg, RhoCOL, Matrix,
                   Pvap25, MaxPvap, Sol25, MW, T25, ChemClass) {
  
  #easy reading kBolts as R
  R = getConst("r")
  
  if (is.na(Kaw25) || Kaw25 == "NA") {
    Kaw25 <- switch (ChemClass,
                     "metal" = 1E-20,
                     "particle" = 1e-20,
                     max(  (ifelse(Pvap25>MaxPvap,MaxPvap,Pvap25) / (Sol25/MW) ) / (R * T25),
                           1E-20) #yet another precaution for too small Kaw and too high Pvap25.
    )
  }
  
  switch(Matrix,
         "air" = 0.54 * (Kow/Kaw25) * Corg * (RhoCOL/1000),
         NA)
  
}

