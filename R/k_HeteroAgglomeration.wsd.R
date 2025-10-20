#' @title Heteroagglomeration (20181102)
#' @name k.HeteroAgglomeration.wsd
#' @description Calculation of the first order rate constant (s-1) for heteroagglomeration of ENPs with other particulates in water and soil
#' largely based on Tufenkji & Elimelech (2003) https://doi.org/10.1021/es034049r
#' @param to.alpha Attachment Efficiency of ENPs with other particulates [-]
#' @param MasConc_Otherparticle Mass concentration of other particulates [kg.m-3]
#' @param from.radius Radius of nanoparticle [m] 
#' @param from.rho Density of nanoparticle [kg.m-3]
#' @param radius_Otherparticle Radius of natural particle [m]
#' @param rho_Otherparticle Density (specific weight) of natural particle [kg/m3]
#' @param rhoFluid Density of fluid matrix [kg/m3]
#' @param Shear Shear rate of the fluid matrix [s-1]
#' @param Temp Temperature of compartment [K]
#' @param DynViscWaterStandard Dynamic viscosity of Water [kg m-1 s-1]
#' @param DynViscAirStandard Dynamic viscosity of Air [kg m-1 s-1]
#' @param SubCompartName Name of relevant subcompartment for which k_HeteroAgglomeration is being calculated
#' @param ScaleName Name of relevant scale for which k_HeteroAgglomeration is being calculated
#' @param COL mass concentration natural colloids in water [kg m-3]
#' @param SUSP mass concentration suspended matter in water [kg m-3]
#' @param RadS Radius of nanoparticle [m]
#' @param RhoS Radius of nanoparticle [kg m-3]
#' @param Shape Shape of the particle [character string]
#' @param rhoMatrix density of the matrix [kg m-3]
#' @param Udarcy Darcy velocity [m s-1]
#' @param hamakerSP.w Hamaker constant for heteroagglomerates [J]
#' @param Matrix type of subcompartment considered 
#' @param RadCOL radius of accumulation mode of particle [m]
#' @param rhoCOL density of accumulation mode of particle [kg m-3]
#' @param rho_CP density of coarse particulate mode of particle [kg m-3]
#' @param radCP radius of accumulation mode of particle [m]
#' @param FRACs fraction of solids in the matrix [-]
#' @param Longest_side description
#' @param Intermediate_side description
#' @param Shortest_side description
#' @param DragMethod Drag method [character string]
#' @return k.HeteroAgglomeration, the rate constant for 1rst order process: heteroagglomeration [s-1]
# #' @seealso \code{\link{f_Brown}}, \code{\link{f_Inter}} and \code{\link{f_Grav}}
#' @export

k_HeteroAgglomeration.wsd <- function(to.alpha,
                                      COL,
                                      SUSP,
                                      Shear,
                                      RadS,
                                      RhoS,
                                      Shape,
                                      RadCOL,
                                      RadCP,
                                      Temp,
                                      DynViscWaterStandard,
                                      RhoCOL,
                                      RhoCP,
                                      rhoMatrix,
                                      Udarcy,
                                      to.FRACs,
                                      hamakerSP.w,
                                      Matrix,
                                      to.SpeciesName,
                                      SubCompartName, 
                                      ScaleName,
                                      Longest_side, Intermediate_side, Shortest_side, 
                                      DragMethod){
  
  if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & 
      (SubCompartName %in% c("agriculturalsoil", "othersoil", "lakesediment", "freshwatersediment"))) {
    return(NA)
  }
  rhoWater = 998 # easyfix could be done more elegantly
  kboltz <- constants::syms$k
  GN <- constants::syms$gn
  
  if(is.na(RadS) || is.null(RadS)) { RadS = Shortest_side/2 }
  
  switch(tolower(Matrix),
          "water" = {
            switch (tolower(to.SpeciesName),
                    "aggregated" = {
                      ColInter <- f_Inter(Shear,RadS,radius_Otherparticle = RadCOL)
                      
                      ColBrown <- f_Brown(Temp=Temp,
                                          viscosity=DynViscWaterStandard,
                                          radius=RadS,
                                          radius_Otherparticle = RadCOL )
                      ColGrav <- f_Grav(rhoParticle= RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        rhoFluid = rhoMatrix,
                                        DynViscWaterStandard = DynViscWaterStandard,
                                        Matrix=Matrix,SubCompartName=SubCompartName, Shape=Shape,
                                        Longest_side=Longest_side, 
                                        Intermediate_side=Intermediate_side, 
                                        Shortest_side=RadS*2, DragMethod=DragMethod)

                      NumConcOther <- f_NumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav+ColInter))
                    },
                    "attached" = {
                      ColInter <- f_Inter(Shear,RadS,radius_Otherparticle = RadCP)
                      
                      ColBrown <- f_Brown(Temp=Temp,
                                          viscosity=DynViscWaterStandard,
                                          radius=RadS,
                                          radius_Otherparticle = RadCP )
                      ColGrav <-f_Grav(rhoParticle= RhoS,
                                       radius_Otherparticle = RadCP,
                                       rho_Otherparticle = RhoCP, 
                                       rhoFluid = rhoMatrix,
                                       DynViscWaterStandard = DynViscWaterStandard,
                                       Matrix=Matrix,SubCompartName=SubCompartName, Shape=Shape,
                                       Longest_side=Longest_side, 
                                       Intermediate_side=Intermediate_side, 
                                       Shortest_side=RadS*2, DragMethod=DragMethod)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCP, 
                                               rho_particle=RhoCP, 
                                               MasConc=SUSP)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav+ColInter))
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
                      ColGrav <- f_Grav(rhoParticle= RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        rhoFluid = rhoWater,
                                        DynViscWaterStandard = DynViscWaterStandard,
                                        Matrix=Matrix,SubCompartName=SubCompartName, Shape=Shape,
                                        Longest_side=Longest_side, 
                                        Intermediate_side=Intermediate_side, 
                                        Shortest_side=RadS*2, DragMethod=DragMethod)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav))
                    },
                    "attached" = {
                      DiffS.w <- f_Diffusivity(Matrix=Matrix, 
                                               Temp = Temp, 
                                               DynVisc=DynViscWaterStandard, 
                                               rad_species=RadS)
                      
                      rhoWater <- 998
                      Por <- 1-to.FRACs
                      GammPDF <- (1-Por)^(1/3)
                      
                      ASPDF <- (2*(1-GammPDF^5))/(2-3*GammPDF+3*GammPDF^5-2*GammPDF^6)
                      aspectratioSFP <- RadS/RadCP
                      PecletNumberFP <- (Udarcy*2*RadCP)/(DiffS.w)
                      vdWaalsNumberSFP <- hamakerSP.w/(kboltz*Temp)
                      
                      BrownSFP <- 2.4*ASPDF^(1/3)*aspectratioSFP^(-0.081)*PecletNumberFP^-0.715*vdWaalsNumberSFP^0.053
                      
                      InterceptSFP <- 0.55*aspectratioSFP^1.55*PecletNumberFP^-0.125*vdWaalsNumberSFP^0.125
                      
                      GravNumberS <- (2*RadS^2*(RhoS-rhoWater)*GN)/(9*DynViscWaterStandard*Udarcy)
                      GravSFP <- 0.22*aspectratioSFP^-0.24*GravNumberS^1.11*vdWaalsNumberSFP^0.053
                      
                      fTotalSFP <- BrownSFP+InterceptSFP+GravSFP
                      
                      Filter <- (3/2)*(1-Por)/(2*RadCP*Por)
                      
                      K_het.sd <- Filter*Udarcy*fTotalSFP*to.alpha 
                      
                      return(K_het.sd)
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
                      ColGrav <- f_Grav(rhoParticle= RhoS,
                                        radius_Otherparticle = RadCOL,
                                        rho_Otherparticle = RhoCOL, 
                                        rhoFluid = rhoWater,
                                        DynViscWaterStandard = DynViscWaterStandard,
                                        Matrix=Matrix,SubCompartName=SubCompartName, Shape=Shape,
                                        Longest_side=Longest_side, 
                                        Intermediate_side=Intermediate_side, 
                                        Shortest_side=RadS*2, DragMethod=DragMethod)
                      
                      NumConcOther <- f_NumConc(rad_particle=RadCOL, 
                                               rho_particle=RhoCOL, 
                                               MasConc=COL)
                      
                      return(to.alpha*NumConcOther*(ColBrown+ColGrav))
                    },
                    "attached" = {
                      DiffS.w <- f_Diffusivity(Matrix=Matrix, 
                                               Temp, DynVisc=DynViscWaterStandard, 
                                               rad_species=RadS)
                      
                      Por <- 1-to.FRACs
                      GammPDF <- (1-Por)^(1/3)
                      
                      ASPDF <- (2*(1-GammPDF^5))/(2-3*GammPDF+3*GammPDF^5-2*GammPDF^6)
                      aspectratioSFP <- RadS/RadCP
                      PecletNumberFP <- (Udarcy*2*RadCP)/(DiffS.w)
                      vdWaalsNumberSFP <- hamakerSP.w/(kboltz*Temp)
                      
                      BrownSFP <- 2.4*ASPDF^(1/3)*aspectratioSFP^(-0.081)*PecletNumberFP^-0.715*vdWaalsNumberSFP^0.053
                      
                      InterceptSFP <- 0.55*aspectratioSFP^1.55*PecletNumberFP^-0.125*vdWaalsNumberSFP^0.125
                      
                      GravNumberS <- (2*RadS^2*(RhoS-rhoWater)*GN)/(9*DynViscWaterStandard*Udarcy)
                      GravSFP <- 0.22*aspectratioSFP^-0.24*GravNumberS^1.11*vdWaalsNumberSFP^0.053
                      
                      fTotalSFP <- BrownSFP+InterceptSFP+GravSFP
                      
                      Filter <- (3/2)*(1-Por)/(2*RadCP*Por)
                      
                      K_het.sd <- Filter*Udarcy*fTotalSFP*to.alpha 
                      
                      return(K_het.sd)
                    },
                    return(NA)
            )
          },
          return(NA)
  )
  

}



