#' @title Area of Marine fraction for each scale
#' @name AreaSea
#' @param all.TotalArea  total AREA of Scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param all.FRACsea  sea fraction of scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param all.ScaleOrder The order from smallest (e.g. local) to largest scale
#' @param all.ScaleNestGroup all scales with the same number are nested within each other
#' @param ScaleName  name of the scale of the box
#' @return Area of Sea for each scale
#' @export
AreaSea <- function (all.TotalArea,
                     all.FRACsea,
                     all.ScaleOrder,
                     all.ScaleNestGroup,
                     ScaleName) {
  
  # Data needed for area calculation
  allData_Scale_Area_Sea <- 
    all.TotalArea |> 
    full_join(all.FRACsea, by = "Scale") |> 
    full_join(all.ScaleOrder, by = "Scale") |> 
    full_join(all.ScaleNestGroup, by = "Scale")
  
  #local function to calculate AreaSea, the naive way (i.e. without "nesting" complications)
  AreaSea4Scale <- function(forScale) {
    AreaSeaData <- NULL
    # ScaleArea <- all.TotalArea$TotalArea[all.TotalArea$Scale == forScale]
    # ScaleFracSea <- all.FRACsea$FRACsea[all.FRACsea$Scale == forScale]
    for(SNG in unique(all.ScaleNestGroup$ScaleNestGroup)){
      AreaSeaData_SNG <-
        allData_Scale_Area_Sea |> 
        filter(ScaleNestGroup == SNG) |> 
        mutate(AreaSea = TotalArea*FRACsea) |> 
        mutate(AreaSea_diff = AreaSea - lead(AreaSea, 
                                             default = 0, 
                                             order_by = -ScaleOrder))
      
      AreaSeaData <- rbind(AreaSeaData,AreaSeaData_SNG)
    }
    
    AreaSeaData |> 
      filter(Scale == forScale) |> 
      pull(AreaSea_diff)
  }
  
  
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
    return(AreaSea4Scale(ScaleName))
  # }
  
}
