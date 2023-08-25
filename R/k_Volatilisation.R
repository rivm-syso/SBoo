#'@title Volatilisation rate constant
#'@name k_Volatilisation
#'@param FRACinf Fraction infiltration #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@returns Volatilisation rate constant [s-1]
#'@export

k_Volatilisation <- function(to.MTC_2w, MTC_2a, to.MTC_2s, Kacompw, FRorig, FRinw, Kscompw,
                             VertDistance, SpeciesName, Matrix, relevant_depth_s, penetration_depth_s){ 
  
  if (SpeciesName %in% c("Molecular")) {
    switch(Matrix,
           "water" = { 
             flux = (to.MTC_2w*MTC_2a/(to.MTC_2w*(Kacompw*FRorig)+MTC_2a))*(Kacompw*FRorig)*FRinw
             return(flux/VertDistance)},
           "soil" = { 
             flux = (to.MTC_2s*MTC_2a)/(to.MTC_2s+MTC_2a/((Kacompw*FRorig)/Kscompw))*f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s)
             return(flux/VertDistance)}
    )
    
    
  } else { 
    return(NA)
  }
}

# #Volatilisation from soil
# #partition coefficient; from excel version:
# #Ksw     (VLOOKUP(ChemClass,QSARtable,7,FALSE))*CORG*RHOsolid/1000
# #Kp.s1R  (FRorig.s1*Ksw+(1-FRorig.s1)*Ksw.alt)*(1000/RHOsolid/CORG)*CORG.s1R
# #Ks1w.R  FRACa.s1R*(Kaw.R*FRorig.s1w)+FRACw.s1R+FRACs.s1R*Kp.s1R*RHOsolid/1000
# #Kaw.R   Kaw*EXP((H0vap/8.314)*(1/298-1/TEMP.R))*EXP(-(H0sol/8.314)*(1/298-1/TEMP.R))*(298/TEMP.R) #but is named KawD.R == v_Kaw()
# 
# #Kp(s,sd) = 
# #Dimensionless soil&sed/water PARTITION COEFFICIENT
# Ksw = FRACa.s1R*(Kaw.R*FRorig.s1w)+FRACw.s1R+FRACs.s1R*Kp.s1R*RHOsolid/1000
# Ksdw = FRACw+FRACs*Kp*RHOsolid/1000
# 
# #from = s or w; to = air
# Volat_s = ((from.MTC_2a*to.MTC_2s)/(from.MTC_2a+to.MTC_2s)/((Kaw*FRorig)/Ksw)) * CORRvolat
# Volat_w = (from.MTC_2a*Kaw.R*FRorig.w1*to.MTC_2w/(from.MTC_2a*Kaw.R*FRorig.w1+to.MTC_2w))*FRw.w0R
# 
# SoilPenetratinDepth = 0.1
# CORRvolat = (EXP((-1/SoilPenetratinDepth)*0)*(1/SoilPenetratinDepth)*DEPTH/(1-EXP((-1/SoilPenetratinDepth)*DEPTH)))
