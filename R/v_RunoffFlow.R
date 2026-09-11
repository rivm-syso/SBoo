#' @title Runoff flow from soil to water
#' @name Runoff
#' @description The runoff is from soil to water. Normally this runs into a river, but, 
#' as there are no rivers modeled in the global scales, it is directed to the sea immediately there. 
#' Basically, it's a flow, but for this reason, modeled as a variable.
#' @param FRACrun volume fraction of precipitation on regional and 
#' continental natural, agricultural and other soil run off to surface water [-] 
#' @param Area area of the subcompartment within the scale
#' @param RAINrate Precipitation for each subcompartment [m/s]
#' @param Compartment Only runoff for Compartment == "soil" 
#' @return Rain water runoff from soils [m3/s]
#' @export
RunoffFlow <- function (FRACrun, Area, RAINrate, Compartment, parent){
  
  out <- Area |>
    dplyr::full_join(RAINrate, by="Scale") |>
    dplyr::full_join(Compartment, by ="SubCompart") |>
    parent$states$clipStates(NoSpeciesKey=TRUE) |>
    dplyr::mutate(
      FRACrun = FRACrun,
      Runoff = dplyr::case_when(
        Compartment == "soil" ~ FRACrun * Area * RAINrate,
        TRUE ~ NA_real_
      )
    ) |>
    dplyr::filter(!is.na(Runoff)) |>
    dplyr::arrange(Scale, SubCompart) |>
    dplyr::select(Scale, SubCompart, Runoff)
  
  return(data.frame(out))
  
  # if (Compartment == "soil") {
  #   return(FRACrun * Area * RAINrate)
  # } else {
  #   return(NA)
  # }
}
