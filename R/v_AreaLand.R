#' @title AreaLand
#' @name AreaLand
#' @description computes the land area in all the scales
#' @param all.TotalArea  total AREA of Scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param all.FRACsea    sea fraction of scale, possibly including other scales,
#' hence the all. It uses nesting to compute the final value.
#' @param all.ScaleOrder The order from smallest (e.g. local) to largest (e.g. Moderate) scale
#' @param all.ScaleNestGroup all scales with the same number are nested within each other
#' @param ScaleName  indicating for which scale the function is called
#' @return Land in SystemArea
#' @export
AreaLand <- function (all.TotalArea,
                      all.FRACsea,
                      all.ScaleOrder,
                      all.ScaleNestGroup,
                      ScaleName) {
  
  # Data needed for area calculation
  allData_Scale_Area_Land <- 
    all.TotalArea |> 
    full_join(all.FRACsea, by = "Scale") |> 
    full_join(all.ScaleOrder, by = "Scale") |> 
    full_join(all.ScaleNestGroup, by = "Scale")
  
  #local function to calculate AreaSea, the naive way (i.e. without "nesting" complications)
  AreaLand4Scale <- function(forScale) {
    AreaLandData <- NULL
    # ScaleArea <- all.TotalArea$TotalArea[all.TotalArea$Scale == forScale]
    # ScaleFracSea <- all.FRACsea$FRACsea[all.FRACsea$Scale == forScale]
    for(SNG in unique(all.ScaleNestGroup$ScaleNestGroup)){
      AreaLandData_SNG <-
        allData_Scale_Area_Land |> 
        filter(ScaleNestGroup == SNG) |> 
        mutate(AreaLand = TotalArea*(1-FRACsea)) |> 
        mutate(AreaLand_diff = AreaLand - lead(AreaLand, 
                                             default = 0, 
                                             order_by = -ScaleOrder))
      
      AreaLandData <- rbind(AreaLandData,AreaLandData_SNG)
    }
    
    AreaLandData |> 
      filter(Scale == forScale) |> 
      pull(AreaLand_diff)
  }
  
    return(AreaLand4Scale(ScaleName))

}
