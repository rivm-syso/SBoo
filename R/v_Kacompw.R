#'@title Dimensionless air/water PARTITION COEFFICIENT of molecular species specific to air compartment (scale)
#'@name Kacompw
#'@description assess the Kaw if it's not given, if Pvap25 is NA or too large, MaxPvap is used
#'@param Pvap25 vapour pressure of original species at 25 C [Pa]
#'@param MaxPvap maximum vapour pressure of original species at 25 C [Pa]
#'@param ChemClass see QSAR table
#'@param Sol25 water solubility of original species at 25 C [mol.m-3]
#'@param H0sol
#'@param T25 const (25 C ==) 298 K
#'@param Tm melting temp
#'@return Kacompw
#'@export
Kacompw <- function(Kaw25, 
                ChemClass, 
                Pvap25, 
                Sol25, 
                H0sol,
                MaxPvap, 
                T25, 
                Tm, Tm_default,
                Temp,
                MW){
  #easy reading kBolts as R
  R = getConst("r")
  #Not in the data
  H0vap = NA
  
  #precaution, cannot be larger than, nor NA
  # CalcPvap25 <- max(ifelse(is.na(Pvap25),MaxPvap,Pvap25), MaxPvap)

  #The actual 
  if (is.na(Kaw25) || Kaw25 == "NA") {
    Kaw25 <- switch (ChemClass,
      "metal" = 1E-20,
      "particle" = 1e-20,
      max(  (ifelse(Pvap25>MaxPvap,MaxPvap,Pvap25) / (Sol25/MW) ) / (R * T25),
          1E-20) #yet another precaution for too small Kaw and too high Pvap25.
    )
  }
  
  if (is.na(H0vap) || H0vap == "NA") { #only used for Kaw
    #this depends on Pvap25, Tm:
    if ((is.na(Tm) || Tm == "NA")) Tm = Tm_default
    H0vap <- 1000 * (-3.82*log(ifelse(Tm>298,Pvap25*exp(-6.79*(1-Tm/T25)),Pvap25))+70)
  }
  
  if (is.na(H0sol) || H0sol == "NA") stop("H0sol is missing")
  if (is.na(Pvap25) || Pvap25 == "NA") stop("Pvap25 is missing")
  Kaw25*exp((H0vap/R)*(1/T25-1/Temp))*exp(-(H0sol/R)*(1/T25-1/Temp))*(T25/Temp)
}


