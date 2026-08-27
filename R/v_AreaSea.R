#' @title AreaSea
#' @name AreaSea
#' @param all.TotalArea  total AREA of Scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param all.FRACsea  sea fraction of scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param ScaleName  name of the scale of the box
#' @return Area of Sea for scale with name ScaleName
#' @export
AreaSea <- function (TotalArea,
                     FRACsea,
                     ScaleName) {
  
  # Nog al vreemde logical in de oude code!
  ContinentalInModerate <- T #For now ! not yet an input option
  
  out <- ScaleName |>
    full_join(TotalArea) |>
    full_join(FRACsea) |>
    mutate(
      AreaSea= TotalArea * (FRACsea)
    ) |>
    mutate(
      AreaSea = dplyr::case_when(
        ScaleName =='Continental' & ContinentalInModerate ~ AreaSea - AreaSea[ScaleName == "Regional"],
        ScaleName == 'Moderate' ~  AreaSea - AreaSea[ScaleName == "Continental"],
        TRUE ~ AreaSea
      )
    ) |>
    select(Scale, AreaSea)
  
  return(out)  
  
  
  # #local function to calculate AreaSea, the naive way (i.e. without "nesting" complications)
  # AreaSea4Scale <- function(forScale) {
  #   ScaleArea <- all.TotalArea$TotalArea[all.TotalArea$Scale == forScale]
  #   ScaleFracSea <- all.FRACsea$FRACsea[all.FRACsea$Scale == forScale]
  #   return(ScaleArea * ScaleFracSea)
  # }
  # 
  # if (ScaleName %in% c("Regional", "Arctic")) {
  #   return(AreaSea4Scale(ScaleName))
  # }
  # 
  # if (ScaleName == "Continental") {
  #   #TotalArea * FRACsea - (TotalArea * FRACsea) for Regional
  #   return(AreaSea4Scale("Continental") - AreaSea4Scale("Regional"))
  # }
  # 
  # ContinentalInModerate <- T #For now ! not yet an input option
  # 
  # if ((ScaleName == "Moderate" & ContinentalInModerate) |
  #     (ScaleName == "Tropic" & !ContinentalInModerate)) {
  #   #TotalArea * FRACsea - (TotalArea * FRACsea) for Continental (including Regional)
  #   return(AreaSea4Scale(ScaleName) - AreaSea4Scale("Continental"))
  # } else {
  #   return(AreaSea4Scale(ScaleName))
  # }
  
}
