#' @title Degradation rate constant measured or calculated
#' @name v_KdegDorC
#' @description calculate k for degradation; 
#' @param FRingas see
#' @param kdeg.air  kdeg.water kdeg.water
#' @param C.OHrad.n ?
#' @param C.OHrad ?
#' @param k0.OHrad ?
#' @param Ea.OHrad ?
#' @param Temp [K]
#' @param T25 = 298 K
#' @return Degradation rate constant for molecular species
#' @export
KdegDorC <- function(kdeg, C.OHrad.n, k0.OHrad, Ea.OHrad, T25, 
                          Q.10, Ksw, Biodeg, CorgStandard, rhoMatrix,
                          Matrix, SpeciesName) {
  if (SpeciesName %in% c("Molecular")) {
    
    switch(Matrix,
           "sediment" = { 
             if(is.na(kdeg)){
               switch (Biodeg,
                       "r" = { a = # following table 5 in report 2015-0161
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,30,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,300,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,3000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,30000,NA))))}, # ready-biodegradable
                       "r-" = { a = 
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,90,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,900,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,9000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
                       "i" = { a = 
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,300000,NA))))}, # inherently biodegradable
                       "p" = { a = 
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,300000,NA)))) }) # persistent
               return(0.1*Q.10^(13/10)*log(2)/a/(3600*24)) } else return(kdeg)
           },
           "soil" = {
             if(is.na(kdeg)){
               switch (Biodeg,
                       "r" = { a = # following table 5 in report 2015-0161
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,30,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,300,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,3000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,30000,NA))))}, # ready-biodegradable
                       "r-" = { a = 
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,90,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,900,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,9000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
                       "i" = { a = 
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,300000,NA))))}, # inherently biodegradable
                       "p" = { a = 
                         ifelse(Ksw/CorgStandard*rhoMatrix/1000<100,300,
                                ifelse(Ksw/CorgStandard*rhoMatrix/1000<1000,3000,
                                       ifelse(Ksw/CorgStandard*rhoMatrix/1000<10000,30000,
                                              ifelse(Ksw/CorgStandard*rhoMatrix/1000>10000,300000,NA)))) }) # persistent
               return(Q.10^(13/10)*log(2)/a/(3600*24)) } else return(kdeg)
           },
           "water" = {
             if(is.na(kdeg)){
               switch (Biodeg,
                       "r" = { Q.10^(13/10)*log(2)/15/(3600*24) }, # ready-biodegradable
                       "r-" = { Q.10^(13/10)*log(2)/50/(3600*24) }, # ready-biodegradable (r-) substances failing the ten-day window
                       "i" = { Q.10^(13/10)*log(2)/150/(3600*24) }, # inherently biodegradable
                       "p" = { 1e-20 }) } else return(kdeg) # persistent
           },
           "air" = {
             if(is.na(kdeg)){
               return(C.OHrad.n * k0.OHrad * exp(-Ea.OHrad/(constants::syms$r*T25)))
             } else return(kdeg)
           })
  } else { 
    switch(Matrix, # particulate
           "air" = kdeg,
           "soil" = kdeg,
           "sediment" = kdeg,
           "water" = kdeg,
           return(NA)
    )
  }
  
}