#' @title Volume of the SubCompartment
#' @name Volume
#' @param VertDistance height or depth for subcompartments below the horizon [m]
#' @param Area area of the compartment [m2]
#' @param FRACcldw fraction of cloudwater [-]
#' @param SubCompartName subcompartment considered
#' @return Volume
#' @export
Volume <- function (VertDistance, Area, FRACcldw, SubCompartName, parent, ScaleName, SpeciesName){
  
  out <- ScaleName |>
    expand_grid(SubCompartName, SpeciesName) |>
    full_join(Area) |>
    full_join(FRACcldw) |>
    full_join(VertDistance, by=c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    mutate(
      FRACcldw = dplyr::case_when(
        SubCompart %in% c('cloudwater', 'air') ~ FRACcldw,
        TRUE ~ NA_real_
      ),
      Volume = dplyr::case_when(
        SubCompart == "air" ~ VertDistance * Area * (1-FRACcldw),
        SubCompart == "cloudwater" ~ VertDistance * Area * FRACcldw,
        TRUE ~ VertDistance * Area
      )
    )  |>
    filter(!is.na(Volume)) |>
    arrange(Scale, SubCompart) |>
    select(Scale, SubCompart, Volume)
  
  return(data.frame(out))
  
  # if(SubCompartName == "air"){
  #   VertDistance * Area * (1-FRACcldw)
  # } else if(SubCompartName == "cloudwater") {
  #   VertDistance * Area * FRACcldw
  # } else 
  #   VertDistance * Area 
}
