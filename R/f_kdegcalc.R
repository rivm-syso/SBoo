#' @title Estimate of degradation rate constant
#' @name kdegcalc
#' @description If the degradation rate in sediment is not mentioned in the substance
#'  data, it can be obtained from (European Commission, 2003a). Following section 3.3.12 of report 2015-0161.
#' @param Ksw 
#' @param Biodeg Biodegradability test result (r / r- / i / p)
#' @param Q.10
#' @param CorgStandard
#' @param RHOsolid
#' @param C.OHrad.n
#' @param k0.OHrad
#' @param Ea.OHrad
#' @return
#' @export
#' 
#'   

kdegcalc <- function (Q.10, Ksw, Biodeg, C.OHrad.n, Ea.OHrad, k0.OHrad,
                    CorgStandard, RHOsolid,
                    Matrix) {
  
  switch(Matrix,
         "sediment" = {
           switch (Biodeg,
                   "r" = { a = # following table 5 in report 2015-0161
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,30,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,300,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,3000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,30000,NA))))}, # ready-biodegradable
                   "r-" = { a = 
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,90,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,900,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,9000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
                   "i" = { a = 
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,300,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,3000,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,30000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,300000,NA))))}, # inherently biodegradable
                   "p" = { a = 
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,300,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,3000,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,30000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,300000,NA)))) }) # persistent
           0.1*Q.10^(13/10)*log(2)/a/(3600*24) },
         "soil" = {
           switch (Biodeg,
                   "r" = { a = # following table 5 in report 2015-0161
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,30,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,300,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,3000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,30000,NA))))}, # ready-biodegradable
                   "r-" = { a = 
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,90,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,900,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,9000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,90000,NA))))}, # ready-biodegradable (r-) substances failing the ten-day window
                   "i" = { a = 
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,300,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,3000,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,30000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,300000,NA))))}, # inherently biodegradable
                   "p" = { a = 
                     ifelse(Ksw/CorgStandard*RHOsolid/1000<100,300,
                            ifelse(Ksw/CorgStandard*RHOsolid/1000<1000,3000,
                                   ifelse(Ksw/CorgStandard*RHOsolid/1000<10000,30000,
                                          ifelse(Ksw/CorgStandard*RHOsolid/1000>10000,300000,NA)))) }) # persistent
           Q.10^(13/10)*log(2)/a/(3600*24) },
         "water" = {
           switch (Biodeg,
                   "r" = { Q.10^(13/10)*log(2)/15/(3600*24) }, # ready-biodegradable
                   "r-" = { Q.10^(13/10)*log(2)/50/(3600*24) }, # ready-biodegradable (r-) substances failing the ten-day window
                   "i" = { Q.10^(13/10)*log(2)/150/(3600*24) }, # inherently biodegradable
                   "p" = { 1e-20 }) # persistent
            },
         "air" = {
           C.OHrad.n * k0.OHrad * exp(-Ea.OHrad/(constants::syms$r*T25))
         }
           )
  

  
}
