#' @title Heteroagglomeration (20181102)
#' @name k.HeteroAgglomeration.wsd
#' @description Calculation of the first order rate constant (s-1) for heteroagglomeration of ENPs with other particulates in water and soil
#' @param alpha Attachment Efficiency of ENPs with other particulates [-]
#' @param MasConc_Otherparticle Mass concentration of other particulates [kg.m-3]
#' @param from.radius Radius of nanoparticle [m] 
#' @param from.rho Density of nanoparticle [kg.m-3]
#' @param radius_Otherparticle Radius of natural particle [m]
#' @param rho_Otherparticle Density (specific weight) of natural particle [kg/m3]
#' @param rhoFluid Density of fluid matrix [kg/m3]
#' @param Shear Shear rate of the fluid matrix [s-1]
#' @param Temp Temperature of compartment [K]
#' @param DynViscWaterStandard Dynamic viscosity of Water []
#' @param DynViscAirStandard Dynamic viscosity of Air []
#' @return k.HeteroAgglomeration, the rate constant for 1rst order process: heteroagglomeration [s-1]
# #' @seealso \code{\link{f_Brown}}, \code{\link{f_Inter}} and \code{\link{f_Grav}}
#' @export
k_HeteroAgglomeration.wsd <- function(alpha,
                                      COL,
                                      SUSP,
                                      Shear,
                                      RadS,
                                      RhoS,
                                      RadCOL,
                                      RadCP,
                                      Temp,
                                      DynViscWaterStandard,
                                      RhoCOL,
                                      RhoCP,
                                      rhoMatrix,
                                      Matrix,
                                      to.SpeciesName,
                                      SubCompartName){
  #for soil and sediment fIntercept assumed 0. Use g in formula, set to 0 in these comaprtments!
  rhoWater = 998 # temp could be done more elegantly
  switch (tolower(Matrix),
          "water" = {
            switch (tolower(to.SpeciesName),
                    "aggregated" = {
                      ColInter <- f_Inter(Shear,RadS,radius_Otherparticle = RadCOL)
                      
                      ColBrown <- f_Brown(Temp=Temp,
                                          viscosity=DynViscWaterStandard,
                                          radius=RadS,
                                          radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(radius = RadS, rho= RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        rhoFluid = rhoMatrix,
                                        DynVisc = DynViscWaterStandard)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(alpha*NumConcOther*(ColBrown+ColGrav+ColInter))
                    },
                    "attached" = {
                      ColInter <- f_Inter(Shear,RadS,radius_Otherparticle = RadCP)
                      
                      ColBrown <- f_Brown(Temp=Temp,
                                          viscosity=DynViscWaterStandard,
                                          radius=RadS,
                                          radius_Otherparticle = RadCP )
                      ColGrav <-f_Grav(radius = RadS, rho= RhoS,
                                       radius_Otherparticle = RadCP,
                                       rho_Otherparticle = RhoCP, 
                                       rhoFluid = rhoMatrix,
                                       DynVisc = DynViscWaterStandard)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCP, 
                                               rho_particle=RhoCP, 
                                               MasConc=SUSP)
                      
                      return(alpha*NumConcOther*(ColBrown+ColGrav+ColInter))
                    },
                    return(NA)
            )
          },
          "soil" = {
            switch (tolower(to.SpeciesName),
                    "aggregated" = {
                      ColBrown <-f_Brown(Temp=Temp,
                                         viscosity=DynViscWaterStandard,
                                         radius=RadS,
                                         radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(radius = RadS, rho= RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        rhoFluid = rhoWater,
                                        DynVisc = DynViscWaterStandard)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(alpha*NumConcOther*(ColBrown+ColGrav))
                    },
                    "attached" = {
                      return(NA) # FP could be integrated here
                    },
                    return(NA)
            )
          },
          "sediment" = {
            switch (tolower(to.SpeciesName),
                    "aggregated" = {
                      ColBrown <- f_Brown(Temp=Temp,
                                          viscosity=DynViscWaterStandard,
                                          radius=RadS,
                                          radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(radius = RadS, rho= RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        rhoFluid = rhoWater,
                                        DynVisc = DynViscWaterStandard)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(alpha*NumConcOther*(ColBrown+ColGrav))
                    },
                    "attached" = {
                      return(NA) # FP could be integrated here
                    },
                    return(NA)
            )
          },
          return(NA)
  )
  
  # ColInter <- f_Inter(Shear,from.radius,radius_Otherparticle )
  # 
  # ColBrown <- f_Brown(Temp,DynVisc,from.radius,radius_Otherparticle )
  # ColGrav <- f_Grav(DynVisc,from.radius,from.rho,radius_Otherparticle ,rho_Otherparticle,g,rhoFluid)
  # 
  # NumConcOther <- f_NumConc(rad_particle=radius_Otherparticle ,rho_particle=rho_Otherparticle, MasConc_Otherparticle)
  # 
  # alpha*NumConcOther*(ColBrown+ColGrav+ColInter)
}



