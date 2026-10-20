#' @title Degradation rate constant measured or calculated
#' @name v_KdegDorC
#' @description calculate k for degradation for particulates and molecules.
#' if no specific value for kdeg is available it is scaled into categories 
#' as defined in Technical Guidance Document on risk assessment in support of
#' Commission Directive 93/67/EEC (European Commission, 2003)
#' @param species mode: which species are needed (Molecular, Small, Large ...) 
#' @param kdeg degradation rate (as input, e.g. based on half life) [s-1]
#' @param C.OHrad OH radical concentration specific to compartment, based on Wania & Daly (2002) [mol m-3]
#' @param C.OHrad.n general OH radical concentration, based on Wania & Daly (2002) [mol m-3]
#' @param k0.OHrad frequency factor of the OH radical reaction [m3 s-1] 
#' @param Ea.OHrad activation energy OH radical reaction [J mol-1]
#' @param T25 298K [K]
#' @param Q.10 rate increase factor per 10C [-]
#' @param KswDorC calculated soil water partitioning coefficient  [-]
#' @param BioDeg biodegradability test result [-]
#' @param CorgStandard Standard mass fraction organic carbon in soil/sediment [-]
#' @param rhoMatrix density of the matrix
#' @param Matrix compartment type considered 
#' @param SpeciesName species name considered
#' @return Degradation rate constant for molecular species
#' @export
KdegDorC <- function(species, Compartment,
                     DegApproach, kdeg, C.OHrad.n, k0.OHrad, Ea.OHrad, T25, 
                     Q.10, KswDorC, Biodeg, CorgStandard, rhoMatrix,
                     Matrix, Shortest_side, Intermediate_side, Longest_side,
                     Kssdr, Shape, CorFacSSA, RadS, degx, degtau, degy, 
                     degtheta, degz, degeta, UVintensity, MICROBconc, Degrading_enzyme) {
  
  # First for Molecular, 
  if ("Unbound" %in% species) {
    Collect_Molecular <- UVintensity |>
      dplyr::left_join(Matrix) |>
      dplyr::left_join(rhoMatrix, by=c("Matrix")) |>
      dplyr::left_join(Compartment) |>
      dplyr::left_join(kdeg) |>
      dplyr::mutate(
        Biodeg = Biodeg,
        tmpBiodeg = KswDorC/CorgStandard*rhoMatrix/1000,
        half_lifes = dplyr::case_when( # Halflife values (table 5 SB 4.0 2015-0161) FOR SOIL SEDIMENT
          Biodeg == 'r' ~ 30, Biodeg == 'i' ~ 90, Biodeg %in% c('i', 'p') ~ 300, TRUE ~ NA),
        half_lifes_corrected_biodegtest = dplyr::case_when(
          tmpBiodeg < 100 ~ half_lifes,  
          tmpBiodeg < 1000 ~ half_lifes * 10,
          tmpBiodeg < 10000 ~ half_lifes * 100,
          tmpBiodeg > 100000 ~ half_lifes * 1000,
          tmpBiodeg < 100000 ~ NA,
          TRUE ~ NA),
        KdegDorC = dplyr::case_when(
          # Soils
          Matrix == 'sediment' & is.na(kdeg) ~ 0.1*Q.10^(13/10)*log(2)/half_lifes_corrected_biodegtest/(3600*24),
          Matrix == 'soil' & is.na(kdeg) ~ Q.10^(13/10)*log(2)/half_lifes_corrected_biodegtest/(3600*24),
          Matrix %in% c('soil','sediment') & !is.na(kdeg) & Biodeg %in% c("r", "r-", "i", "p") ~ kdeg,
          # Water
          Matrix == 'water' & is.na(kdeg) & Biodeg == 'r' ~ Q.10^(13/10)*log(2)/15/(3600*24),
          Matrix == 'water' & is.na(kdeg) & Biodeg == 'r-' ~ Q.10^(13/10)*log(2)/50/(3600*24),
          Matrix == 'water' & is.na(kdeg) & Biodeg == 'i' ~ Q.10^(13/10)*log(2)/150/(3600*24),
          Matrix == 'water' & is.na(kdeg) & Biodeg == 'p' ~ 1e-20,
          Matrix == 'water' & !is.na(kdeg) & Biodeg %in% c("r", "r-", "i", "p") ~ kdeg,
          # air
          Matrix == 'air' & is.na(kdeg) ~ C.OHrad.n * k0.OHrad * exp(-Ea.OHrad/(constants::syms$r*T25)),
          Matrix == 'air' & !is.na(kdeg) ~ kdeg
        ),
        Species = "Unbound"
      ) |>
      dplyr::select(Scale,SubCompart,Species,KdegDorC)
  }
  
  # For particulate
  if (any(c("Large", "Small", "Solid") %in% species)) {
    
    #Set default shortest, intermediate and longest side, in case it is not defined
    if (is.na(Shortest_side) || is.null(Shortest_side) ) {
      Shortest_side <- RadS * 2
    }
    if ( is.na(Intermediate_side) || is.null(Intermediate_side) ) {
      Intermediate_side <- RadS * 2
    }
    
    if ( is.na(Longest_side) || is.null(Longest_side) ) {
      Longest_side <- RadS * 2
    }
    if (!DegApproach %in% c("Default","Kssdr","PlasticFADE")){
      warning("v_KdegDorC: DegApproach not set to valid option, using Default.")
      DegApproach <- "Default"
    }
    if(is.null(Degrading_enzyme) || is.na(Degrading_enzyme)){
      warning("v_KdegDorC: Degrading_enzyme is not set (NA), now set to FALSE")
      Degrading_enzyme = "FALSE"
    }
  
    Collect_Particulate <- UVintensity |>
      dplyr::left_join(Matrix) |>
      dplyr::left_join(MICROBconc) |>
      dplyr::left_join(Compartment) |>
      dplyr::left_join(kdeg) |>
      dplyr::mutate(
        SAV = dplyr::case_when(
          !is.na(Shape) & Shape %in% c("Cylindric - circular", "Fiber") ~ 4/(Shortest_side*100)+2/(Longest_side*100),
          !is.na(Shape) & Shape %in% c("Film") ~ 2/(Shortest_side*100)+2/(Intermediate_side*100)+2/(Longest_side*100),
          TRUE ~ 6/(Shortest_side*100)),
        ShapeFac = dplyr::case_when(
          !is.na(Shape) & Shape %in% c("Cylindric - circular", "Fiber") ~ 3,
          !is.na(Shape) & Shape %in% c("Film") ~ 2,
          TRUE ~ 4),
        UVintensity = dplyr::case_when(
          is.na(UVintensity) & DegApproach == "PlasticFADE" ~ 0,
          TRUE ~ UVintensity),
        Degrading_enzyme = dplyr::case_when(
          DegApproach == 'PlasticFADE' & is.na(Degrading_enzyme) ~ "FALSE",
          TRUE ~ Degrading_enzyme),
        properMatrix = Matrix %in% c("air", "soil", "sediment", "water")
      ) 
    
    Collect_Particulate <- Collect_Particulate |>
      dplyr::mutate(
        KdegDorC = dplyr::case_when(
          DegApproach %in% c("Kssdr", "PlasticFADE", "Default") & 
            (!is.na(Degrading_enzyme) & (Degrading_enzyme == "FALSE" | Degrading_enzyme == FALSE) & (!is.na(UVintensity) & UVintensity == 0)) ~
            0,
          DegApproach == "PlasticFADE" ~ degx * SAV^degtau * (degy * UVintensity^degtheta + degz * MICROBconc^degeta) / (24*60*60),
          DegApproach == "Kssdr" ~ ShapeFac * Kssdr * CorFacSSA / (Shortest_side / 2),
          DegApproach == "Default" & properMatrix ~ kdeg,
          # ? DegApproach == "Default" & !Matrix %in% c("air", "soil", "sediment", "water") ~ NA,
          TRUE ~ NA
        )
      ) |> 
      dplyr::filter(!is.na(KdegDorC)) |>
      dplyr::select(Scale,SubCompart,KdegDorC) |>
      # for all Species except Unbound
      tidyr::crossing(Species = species[species != "Unbound"])
  }
  
  out <-
    if (!is.null(Collect_Molecular) && !is.null(Collect_Particulate)) {
      dplyr::bind_rows(Collect_Molecular, Collect_Particulate)
    } else if (!is.null(Collect_Molecular)) {
      Collect_Molecular
    } else if (!is.null(Collect_Particulate)) {
      Collect_Particulate
    } else {
      NULL
    }
  
  return(data.frame(out))
}
  
  # #Set default shortest, intermediate and longest side, in case it is not defined
  # if ( is.na(Shortest_side) || is.null(Shortest_side) ) {
  #   Shortest_side <- RadS * 2
  # }
  # 
  # if ( is.na(Intermediate_side) || is.null(Intermediate_side) ) {
  #   Intermediate_side <- RadS * 2
  # }
  # 
  # if ( is.na(Longest_side) || is.null(Longest_side) ) {
  #   Longest_side <- RadS * 2
  # }
  # # browser()
  # if (!DegApproach %in% c("Default","Kssdr","PlasticFADE")){
  #   warning("v_KdegDorC: DegApproach not set to valid option, using Default.")
  #   DegApproach <- "Default"
  # }
  # 
  # if (SpeciesName %in% c("Molecular")) {
  #   
  #   switch(Matrix,
  #          "sediment" = { 
  #            if(is.na(kdeg)){
  #              switch (Biodeg,
  #                      "r" = { a = # following table 5 in report 2015-0161
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,30,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,300,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,3000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,30000,NA))))}, # ready-biodegradable
  #                      "r-" = { a = 
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,90,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,900,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,9000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
  #                      "i" = { a = 
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA))))}, # inherently biodegradable
  #                      "p" = { a = 
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA)))) }) # persistent
  #              return(0.1*Q.10^(13/10)*log(2)/a/(3600*24)) } else return(kdeg)
  #          },
  #          "soil" = {
  #            if(is.na(kdeg)){
  #              switch (Biodeg,
  #                      "r" = { a = # following table 5 in report 2015-0161
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,30,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,300,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,3000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,30000,NA))))}, # ready-biodegradable
  #                      "r-" = { a = 
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,90,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,900,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,9000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
  #                      "i" = { a = 
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA))))}, # inherently biodegradable
  #                      "p" = { a = 
  #                        ifelse(KswDorC/CorgStandard*rhoMatrix/1000<100,300,
  #                               ifelse(KswDorC/CorgStandard*rhoMatrix/1000<1000,3000,
  #                                      ifelse(KswDorC/CorgStandard*rhoMatrix/1000<10000,30000,
  #                                             ifelse(KswDorC/CorgStandard*rhoMatrix/1000>100000,300000,NA)))) }) # persistent
  #              return(Q.10^(13/10)*log(2)/a/(3600*24)) } else return(kdeg)
  #          },
  #          "water" = {
  #            if(is.na(kdeg)){
  #              switch (Biodeg,
  #                      "r" = { Q.10^(13/10)*log(2)/15/(3600*24) }, # ready-biodegradable
  #                      "r-" = { Q.10^(13/10)*log(2)/50/(3600*24) }, # ready-biodegradable (r-) substances failing the ten-day window
  #                      "i" = { Q.10^(13/10)*log(2)/150/(3600*24) }, # inherently biodegradable
  #                      "p" = { 1e-20 }) } else return(kdeg) # persistent
  #          },
  #          "air" = {
  #            if(is.na(kdeg)){
  #              return(C.OHrad.n * k0.OHrad * exp(-Ea.OHrad/(constants::syms$r*T25)))
  #            } else return(kdeg)
  #          })
  # } else { 
  #   
  #   #Determine the shape factor, based on the approach of Maga et al. (2022) for surface degradation rate
  #   if (!is.na(Shape) && (Shape == "Cylindric - circular" | Shape == "Fiber")) {
  #     ShapeFac <- 3
  #     SAV <- 4/(Shortest_side*100)+2/(Longest_side*100)#Surface area to volume ratio - in cm-1
  #   } else if (!is.na(Shape) && Shape == "Film") {
  #     ShapeFac <- 2
  #     SAV <- 2/(Shortest_side*100)+2/(Intermediate_side*100)+2/(Longest_side*100)
  #   } else {
  #     ShapeFac <- 4 #If the shape is sphere, not given or irregular, use the Kssdr calculation for a sphere
  #     SAV <- 6/(Shortest_side*100)
  #   }
  #   
  #   switch(DegApproach,     # Calculate kdeg (s-1) using either of the following 3 approaches
  #          "PlasticFADE" = {
  #            
  #            if(is.na(UVintensity)){
  #              warning("v_KdegDorC: UVintensity is not set (NA), now set to 0")
  #              UVintensity = 0
  #            }
  #            if(is.na(Degrading_enzyme)){
  #              warning("v_KdegDorC: Degrading_enzyme is not set (NA), now set to FALSE")
  #              Degrading_enzyme = "FALSE"
  #            }
  #            
  #            #Some polymers cannot be degraded in the absence of UV light (e.g. polyolefins) (UVintensity=0), if degrading enzymes don't exist in the environment for the specific polymer chain (Degrading_enzyme=FALSE). 
  #            #kdeg evaluates to 0 in that compartment for that polymer in that case.
  #            #However, these polymers can be degraded in the presence of UV, as UV initiate the breakdown of the polymer chainl, allowing other enzymes to degrade the polymer.
  #            
  #            if (!is.na(Degrading_enzyme) &&
  #                (Degrading_enzyme == "FALSE" || Degrading_enzyme == FALSE) &&
  #                (!is.na(UVintensity) && UVintensity == 0)) {
  #               kdeg <- 0
  #             } else {
  #            kdeg <- degx * SAV^degtau * (degy * UVintensity^degtheta + degz * MICROBconc^degeta) / (24*60*60)
  #            }
  #          },
  #          
  #          "Kssdr" = {
  #            #With the SSDR approach, k deg=0 when UV=0 and Degrading enzymes don't exist as well (same as for plasticFADE model)
  #            if (!is.na(Degrading_enzyme) &&
  #                (Degrading_enzyme == "FALSE" || Degrading_enzyme == FALSE) &&
  #                (!is.na(UVintensity) && UVintensity == 0)) {
  #              kdeg <- 0
  #            } else {
  #              kdeg <- ShapeFac * Kssdr * CorFacSSA / (Shortest_side / 2)
  #            }
  #          },
  #          
  #          "Default" = {
  #            switch(Matrix, #particulate
  #                   "air" = kdeg,
  #                   "soil" = kdeg,
  #                   "sediment" = kdeg,
  #                   "water" = kdeg,
  #                   NA)
  #          },
  #          NA
  #   )
  #   
  # }
  # 
