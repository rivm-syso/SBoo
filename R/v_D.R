#'@title D
#'@name D
#'@description Apparent octanol/water PARTITION COEFFICIENT at neutral pH
#'@param ChemClass
#'@param FRorig 
#'@param pKa
#'@param Kow
#'@return 
#'@export
D <- function(FRorig, pKa, Kow, ChemClass){
  
  Kow.alt = 10^(log10(Kow)-3.5)
  
  
  switch(ChemClass,
         "acid" = 
           1/(1+10^(7-pKa)) * Kow + (1 - 1/(1+10^(7-pKa)))*Kow.alt,
         "base" = 
           1/(1+10^(pKa-7)) * Kow + (1- 1/(1+10^(pKa-7)))*Kow.alt,
         # else for other ChemClass:
         Kow
          
       )
}
