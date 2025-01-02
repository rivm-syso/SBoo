
f_Soil.wetweight <- function(Conc.soil, # in kg/m3 soil or sediment
                             Fracw, Fracs, RHOsolid, RHOw){
  Conc.soil*1000/(Fracw*RHOw+Fracs*RHOsolid) # in g/kg (wet) soil
}
