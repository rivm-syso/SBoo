#' @title Rate constant for wet deposition of particulate species due to rain
#' @name k_WetDeposition
#' @description Calculation of the first order rate constant for particle deposition from cloudwater compartment to the soil or water surface [s-1]
#' @param to.AREA Surface area of receiving land or water compartment [m2]
#' @param from.VOLUME Volume of cloudwater compartment [m3]
#' @param RAINrate average rain rate [m s-1]
#' @export
k_WetDeposition <- function (Area, Volume, RAINrate, ScaleName, SubCompartName, SpeciesName, parent){
  
  out <- parent$FromDataAndTo("k_WetDeposition")
  # Names of dataframe
  names <- dplyr::all_of(names(out))

  # Build
  VolumeCloudwater <- Volume |>
    dplyr::filter(SubCompart == "cloudwater") |>
    dplyr::select(Scale, Volume)
  out <- out |>
    dplyr::left_join(VolumeCloudwater) |>
    dplyr::left_join(RAINrate) |>
    # Area of the to-Subcompart 
    dplyr::left_join(Area, by = c("to.SubCompart" = "SubCompart", "Scale")) |>
    dplyr::filter(!is.na(Area)) |>

  # Calc

    dplyr::mutate(
      k_WetDeposition = (RAINrate* Area)/Volume
    ) |>
    dplyr::filter(!is.na(k_WetDeposition)) |>
    dplyr::select(names, k_WetDeposition) |>
    dplyr::arrange(dplyr::across(names))

  return(data.frame(out))
  
  # (RAINrate* to.Area)/from.Volume
}

