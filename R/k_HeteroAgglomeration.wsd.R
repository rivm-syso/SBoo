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
#' @param DynVisc Dyanic DynVisc of fluid matrix []
#' @return k.HeteroAgglomeration, the rate constant for 1rst order process: heteroagglomeration [s-1]
# #' @seealso \code{\link{f_Brown}}, \code{\link{f_Inter}} and \code{\link{f_Grav}}
#' @export
k_HeteroAgglomeration.wsd <- function(to.alpha,
                                      COL,
                                      SUSP,
                                      Shear,
                                      from.RadS,
                                      from.RhoS,
                                      RadCOL,
                                      RadCP,
                                      Temp,DynVisc,
                                      RhoCOL,
                                      RhoCP,
                                      rhoMatrix,
                                      Matrix,
                                      to.Species){
  #for soil and sediment fIntercept assumed 0. Use g in formula, set to 0 in these comaprtments!
  
  switch (Matrix,
          "Water" = {
            switch (to.Species,
                    "Small" = {
                      ColInter <- f_Inter(Shear,from.RadS,radius_Otherparticle = RadCOL)
                      
                      ColBrown <- f_Brown(Temp,DynVisc,from.RadS,radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(DynVisc,from.RadS, from.RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        g, rhoMatrix)
                      
                      NumConcOther <- fNumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav+ColInter))
                    },
                    "Large" = {
                      ColInter <- f_Inter(Shear,from.RadS,radius_Otherparticle = RadCP)
                      
                      ColBrown <- f_Brown(Temp,DynVisc,from.RadS,radius_Otherparticle = RadCP )
                      ColGrav <- f_Grav(DynVisc,from.RadS, from.RhoS,
                                        radius_Otherparticle = RadCP,
                                        rho_Otherparticle = RhoCP, 
                                        g, rhoMatrix)
                      
                      NumConcOther <- fNumConc(rad_particle=RadCP, 
                                               rho_particle=RhoCP, 
                                               MasConc=SUSP)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav+ColInter))
                    },
                    return(NA)
            )
          },
          "Soil" = {
            switch (to.Species,
                    "Small" = {
                      ColBrown <- f_Brown(Temp,DynVisc,from.RadS,radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(DynVisc,from.RadS, from.RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        g, rhoMatrix)
                      
                      NumConcOther <- fNumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav))
                    },
                    "Large" = {
                      return(NA) # FP could be integrated here
                    },
                    return(NA)
            )
          },
          "Sediment" = {
            switch (to.Species,
                    "Small" = {
                      ColBrown <- f_Brown(Temp,DynVisc,from.RadS,radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(DynVisc,from.RadS, from.RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        g, rhoMatrix)
                      
                      NumConcOther <- fNumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav))
                    },
                    "Large" = {
                      return(NA) # FP could be integrated here
                    },
                    return(NA)
            )
          },
          return(NA)
  )
  
  ColInter <- f_Inter(Shear,from.radius,radius_Otherparticle )
  
  ColBrown <- f_Brown(Temp,DynVisc,from.radius,radius_Otherparticle )
  ColGrav <- f_Grav(DynVisc,from.radius,from.rho,radius_Otherparticle ,rho_Otherparticle,g,rho.fluid)
  
  NumConcOther <- fNumConc(rad_particle=radius_Otherparticle ,rho_particle=rho_Otherparticle, MasConc_Otherparticle)
  
  to.alpha*NumConcOther*(ColBrown+ColGrav+ColInter)
}



