#' @title AreaLand
#' @name AreaLand
#' @description computes the land area in all the scales
#' @param all.TotalArea  total AREA of Scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param all.FRACsea    sea fraction of scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param ScaleName  indicating for which scale the function is called
#' @return Land in SystemArea
#' @export
AreaLand <- function (TotalArea,
                      FRACsea,
                      ScaleName) {
  
  # Nog al vreemde logical in de oude code!
  ContinentalInModerate <- T #For now ! not yet an input option
  
  out <- ScaleName |>
    full_join(TotalArea) |>
    full_join(FRACsea) |>
    mutate(
      AreaLand = TotalArea * (1-FRACsea)
    ) |>
    mutate(
      AreaLand = dplyr::case_when(
        ScaleName =='Continental' & ContinentalInModerate ~ AreaLand - AreaLand[ScaleName == "Regional"],
        ScaleName == 'Moderate' ~  AreaLand - AreaLand[ScaleName == "Continental"],
        TRUE ~ AreaLand
      )
    ) |>
    select(Scale, AreaLand)
  
  return(out)  
  
  # #local function to calculate AreaSea, the naive way (i.e. without "nesting" complications)
  # AreaLand4Scale <- function(forScale) {
  #   ScaleArea <- all.TotalArea$TotalArea[all.TotalArea$Scale == forScale]
  #   ScaleFracLand <- 1 - all.FRACsea$FRACsea[all.FRACsea$Scale == forScale]
  #   return(ScaleArea * ScaleFracLand)
  # }
  # 
  # if (ScaleName %in% c("Regional", "Arctic")) {
  #   return(AreaLand4Scale(ScaleName))
  # }
  # 
  # if (ScaleName == "Continental") {
  #   #TotalArea * FRACLand - (TotalArea * FRACLand) for Regional
  #   return(AreaLand4Scale("Continental") - AreaLand4Scale("Regional"))
  # }
  # 
  # ContinentalInModerate <- T #For now ! not yet an input option
  # 
  # if ((ScaleName == "Moderate" & ContinentalInModerate) |
  #     (ScaleName == "Tropic" & !ContinentalInModerate)) {
  #   #TotalArea * FRACLand - (TotalArea * FRACLand) for Continental (including Regional)
  #   return(AreaLand4Scale(ScaleName) - AreaLand4Scale("Continental"))
  # } else {
  #   return(AreaLand4Scale(ScaleName))
  # }

}
