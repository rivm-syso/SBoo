#' @title Fragmentation of macroplastics to microplastics
#' @name k_MacroFragmentation
#' @description To Be Completed
#' @param kMfrag Falling appart or fragmenting rate constant [s-1]
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_MacroFragmentation [s-1]
#' @export


KfragDorC <- function (kMfrag, SubCompartName, Matrix,
                       FragApproach,
                       ScaleName,
                       Regional_and_Continental_deepocean,
                       Shortest_side, Intermediate_side, Longest_side,
                       WINDspeed, Volume, RhoS,
                       Kssdr, Shape, CorFacSSA, RadS, degx, degtau, degy, 
                       UVintensity, VeloMedium, CORdrag,
                       CoefDrag,
                       coeff_of_friction,
                       power_type,
                       fraga, fragdelta, fragb, fragalpha, fragc, fragabeta
) {
  
  if (((ScaleName %in% c("Tropic", "Moderate", "Arctic")) &
       (
         SubCompartName == "freshwatersediment" |
         SubCompartName == "lakesediment" |
         SubCompartName == "lake" |
         SubCompartName == "river" |
         SubCompartName == "agriculturalsoil" |
         SubCompartName == "othersoil"
       )
  )) {
    return(NA)
  } else if (!is.null(ScaleName) &&
             ScaleName %in% c("Regional", "Continental") &&
             (SubCompartName == "deepocean") &&
             (
               is.na(Regional_and_Continental_deepocean) ||
               isFALSE(Regional_and_Continental_deepocean) ||
               Regional_and_Continental_deepocean == "FALSE"
             )) {
    return(NA)
  } else {
    #Set default shortest, intermediate and longest side, in case it is not defined
    if ( is.na(Shortest_side) || is.null(Shortest_side) ) {
      Shortest_side <- RadS * 2
    }
    
    if ( is.na(Intermediate_side) || is.null(Intermediate_side) ) {
      Intermediate_side <- RadS * 2
    }
    
    if ( is.na(Longest_side) || is.null(Longest_side) ) {
      Longest_side <- RadS * 2
    }
    #Determine the shape factor, based on the approach of Maga et al. (2022) for surface degradation rate
    if (!is.na(Shape) && (Shape == "Cylindric - circular" | Shape == "Fiber")) {
      ShapeFac <- 3
      SAV <- 4/(Shortest_side*100)+2/(Longest_side*100)#Surface area to volume ratio - in cm-1
      ContArea = Shortest_side * Longest_side # ContArea: Contact Area 
    } else if (!is.na(Shape) && Shape == "Film") {
      ShapeFac <- 2
      SAV <- 2/(Shortest_side*100)+2/(Intermediate_side*100)+2/(Longest_side*100)
      ContArea = Intermediate_side * Longest_side # ContArea: Contact Area 
    } else {
      ShapeFac <- 4 #If the shape is sphere, not given or irregular, use the Kssdr calculation for a sphere
      SAV <- 6/(Shortest_side*100)
      ContArea = 2 * pi * Shortest_side # ContArea: Contact Area 
    }
    
    switch (SubCompartName,
            "air" = {vMedium = WINDspeed},  # m/s
            {vMedium = VeloMedium}
    )
    
    #### FragApproach ####
    switch(FragApproach,     # Calculate kdeg (s-1) using either of the following 3 approaches
           "PlasticFADE" = {
             # if(SubCompartName == "river") {browser()}
             # browser()
             if(is.na(UVintensity)){
               warning("v_KdegDorC: UVintensity is not set (NA), now set to 0")
               UVintensity = 0
             }
             
             switch(power_type,
                    "None" = {MechPower = 0},
                    "Drag" = {
                      
                      #  CORdrag is fc (f_C) which is size dependent correction factor for estimatin the plastic object velocity due to drag.
                      MechPower = 0.5 * 1000 * CoefDrag*ContArea*(vMedium*CORdrag)^3 # in mW
                    },
                    "Friction" = {
                      GN <- constants::syms$gn
                      ObjectVolume = fVol(rad_particle=rad_particle,Shape=Shape, Longest_side=Longest_side, Intermediate_side=Intermediate_side, Shortest_side=Shortest_side)
                      ObjMass_kg = ObjectVolume * RhoS
                      MechPower = coeff_of_friction * ObjMass_kg * GN * vMedium * CORdrag
                    },
                    "NA" = {MechPower = 0},
                    NA
             )
             kMfragCalc <- fraga * SAV^fragdelta * (fragb * UVintensity^fragalpha + fragc * MechPower^fragabeta)/ (24*60*60)
             return(kMfragCalc)
             # }
           },
           
           "Kssdr" = {
             #With the SSDR approach, k deg=0 when UV=0 and Degrading enzymes don't exist as well (same as for plasticFADE model)
             if ((!is.na(UVintensity) && UVintensity == 0)) {
               kMfragCalc <- 0
             } else {
               kMfragCalc <- ShapeFac * Kssdr * CorFacSSA / (Shortest_side / 2)
               return(kMfragCalc)
             }
           },
           "Default" = {
             switch(Matrix, #particulate
                    "air" =   return(kMfrag),
                    "soil" =   return(kMfrag),
                    "sediment" =   return(kMfrag),
                    "water" =   return(kMfrag),
                    NA)
           },
           NA)
  }
  

  
}