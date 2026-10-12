#' @title Rate constant for wet deposition of particulate species due to rain
#' @name k_WetDeposition
#' @description Calculation of the first order rate constant for particle deposition from cloudwater compartment to the soil or water surface [s-1]
#' @param to.AREA Surface area of receiving land or water compartment [m2]
#' @param from.VOLUME Volume of cloudwater compartment [m3]
#' @param RAINrate average rain rate [m s-1]
#' @export
k_WetDeposition <- function (Area, Volume, RAINrate, ScaleName, SubCompartName, SpeciesName, parent){
  
  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_WetDeposition")))

  # Build grid
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Volume", "RAINrate")
  dfs <- mget(args, envir = environment())
  grid_from <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  args <- c("ScaleName", "SubCompartName", "SpeciesName", "Area")
  dfs <- mget(args, envir = environment())
  grid_to <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()

  # Output
  out <- parent$FromDataAndTo("k_WetDeposition") |>
    # From
    dplyr::left_join(grid_from, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species")) |>

    #To
    dplyr::left_join(grid_to, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species"), suffix=c("_from", "_to")) |>

    dplyr::mutate(
      k_WetDeposition = (RAINrate* Area)/Volume
    ) |>
    dplyr::filter(!is.na(k_WetDeposition)) |>
    dplyr::select(names, k_WetDeposition) |>
    dplyr::arrange(dplyr::across(names))

  return(data.frame(out))
  
  # (RAINrate* to.Area)/from.Volume
}

