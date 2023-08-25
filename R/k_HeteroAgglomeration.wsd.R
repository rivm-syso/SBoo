#' @title Heteroagglomeration (20181102)
#' @name k.HeteroAgglomeration.wsd
#' @description Calculation of the first order rate constant (s-1) for heteroagglomeration of ENPs with other particulates in water and soil
#' @param to.alpha Attachment Efficiency of ENPs with other particulates [-]
#' @param MasConc_Otherparticle Mass concentration of other particulates [kg.m-3]
#' @param from.radius Radius of nanoparticle [m] 
#' @param from.rho Density of nanoparticle [kg.m-3]
#' @param radius_Otherparticle Radius of natural particle [m]
#' @param rho_Otherparticle Density (specific weight) of natural particle [kg/m3]
#' @param rho.fluid Density of fluid matrix [kg/m3]
#' @param Shear Shear rate of the fluid matrix [s-1]
#' @param Temp Temperature of compartment [K]
#' @param viscosity Dyanic viscosity of fluid matrix []
#' @return k.HeteroAgglomeration, the rate constant for 1rst order process: heteroagglomeration [s-1]
#' @seealso \code{\link{fBrown}}, \code{\link{fInter}} and \code{\link{fGrav}}
#' @export
k.HeteroAgglomeration.wsd <- function(to.alpha,
                                      MasConc_Otherparticle,
                                      Shear,
                                      from.radius,
                                      radius_Otherparticle ,
                                      Temp,viscosity,
                                      from.rho,
                                      rho_Otherparticle,
                                      rho.fluid){
  #for soil and sediment fIntercept assumed 0. Use g in formula, set to 0 in these comaprtments!
  ColInter <- fInter(Shear,from.radius,radius_Otherparticle )
  ColBrown <- fBrown(Temp,viscosity,from.radius,radius_Otherparticle )
  ColGrav <- fGrav(viscosity,from.radius,from.rho,radius_Otherparticle ,rho_Otherparticle,g,rho.fluid)
  
  NumConcOther <- fNumConc(rad_particle=radius_Otherparticle ,rho_particle=rho_Otherparticle, MasConc_Otherparticle)
  
  to.alpha*NumConcOther*(ColBrown+ColGrav+ColInter)
}



