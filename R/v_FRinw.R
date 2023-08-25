#'@title FRaction of chemical truely dissolved in water phase fresh water (relevant to Molecular species)
#'@name Frinw
#'@description not the part suspended or colloidal 
#'@param Kp.susp
#'@param D , to calculate Kp.col
#'@param SUSP
#'@param SUSP
#'@return 
#'@export
FRinw <- function(FRorig_spw, FRACw, FRACa, FRACs, Kp, RHOsolid, KpCOL, 
                  Kacompw, SUSP, COL, Matrix){
  #TODO next formula in f_ function
  
  switch(Matrix,
         "water" = 1/(1+Kp*SUSP/1000+KpCOL*COL/1000),
         "soil" = # for soil pore water
           FRACw/(FRACa*(Kacompw*FRorig_spw)+FRACw+FRACs*Kp*RHOsolid/1000)
  )
}
