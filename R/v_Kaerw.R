#' @title Aerosol water - air Partition coefficient
#' @name Kaerw
#' @description pm
#' @param Kacompw  air water partitioning coefficient [-]
#' @param FRorig fraction of original species [-]
#' @param SubCompartName subcompartment considered
#' @return Kaerw
#' @export
Kaerw <- function (Kacompw, FRorig, SubCompartName) {
  
  out <- Kacompw |>
    expand_grid(SubCompartName) |>
    full_join(FRorig, by="SubCompart") |>
    mutate(
      Kaerw = dplyr::case_when(
        SubCompartName == 'air' ~ 1/(Kacompw*FRorig),
        TRUE ~ NA_real_
      )
    ) |>
    filter(!is.na(Kaerw)) |>
    arrange(Scale, SubCompart) |>
    select(Scale, SubCompart, Kaerw)
  
  return(out)
   
  # switch(SubCompartName,
  #        "air" = 1/(Kacompw*FRorig),
  #        NA)
}