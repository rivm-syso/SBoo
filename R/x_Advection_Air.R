#' @title Advection of fluxes of air between scales 
#' @name x_Advection_air
#' @description Calculation using AirFlow, but adjust to mass balance volume flows
#' @param AirFlow to pick for the "from" scale 
#' @param ScaleName name of the scale of the source/sink of the process
#' @param Remove_global If this variable is TRUE, the global scales (Arctic, Moderate and Tropic) are removed
#' @return x.Advection_Air, the flow
#' @export
x_Advection_Air <- function(AirFlow, Remove_global, parent){
  # The airflow was step 1; now the corrections for the mass balance to add up..
  # "normal" k is 1/Tau ;  Airflow = Volume/Tau
  
  out <- parent$FromDataAndTo("x_Advection_Air")
  
  # to keep mass balance
  Cont2Regional <- AirFlow$AirFlow[AirFlow$Scale == "Regional"]
  Moder2Cont <- AirFlow$AirFlow[AirFlow$Scale == "Continental"] - Cont2Regional
  
  out <- out |>
    dplyr::left_join(AirFlow, by = c("from.Scale" = "Scale", "SubCompart")) |>
    dplyr::rename(from.AirFlow = AirFlow) |>
    dplyr::left_join(AirFlow, by = c("to.Scale" = "Scale", "SubCompart")) |>
    dplyr::rename(to.AirFlow = AirFlow) |>
    
    dplyr::mutate(flow = dplyr::case_when(
      #exceptions first
      from.Scale == "Continental" & to.Scale == "Regional"   ~ Cont2Regional,  # scalar or length n
      from.Scale == "Moderate"    & to.Scale == "Continental" ~ Moder2Cont,
      from.Scale == "Continental" & to.Scale == "Moderate"    ~ Moder2Cont,
      from.Scale %in% c("Arctic", "Regional", "Tropic")       ~ from.AirFlow,
      TRUE                                                     ~ to.AirFlow
    )) |> 
    
    # remove columns
    select(-from.AirFlow, -to.AirFlow)
  
  return (out)
}