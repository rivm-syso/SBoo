#'@title Desorption of molecular species from sediment to water
#'@name k_Desorption 
#'@description desorption based on the two-film resistance model based on Schwarzenbach et al. (1993) ISBN: 978-1-118-76723-8
#'@param Ksdcompw sediment water partitioning coefficient, see Ksdcompw [-]
#'@param MTC_2w partial mass transfer coefficient to water, see MTC_2w [m s-1]
#'@param MTC_2sd partial mass transfer coefficient to sediment, see MTC_2sd [m s-1]
#'@param SpeciesName name of the species considered
#'@param SubCompartName subcompartment considered
#'@param ScaleName scale considered
#'@param Test Test = TRUE mimics SB4 in Excel version, Test = FALSE includes SB enhancements
#'@param VertDistance vertical distance of compartment [m]
#'@return Desorption rate constant from sediment to water [s-1]
#'@export
k_Desorption <- function (Ksdcompw, MTC_2w, MTC_2sd, VertDistance,
                          SpeciesName, to.SubCompartName, SubCompartName, ScaleName, Test, parent) {
  
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Desorption")))
  
  # Build grid
  args <- c("SubCompartName", "ScaleName", "SpeciesName", "VertDistance", "Ksdcompw", "MTC_2w")
  dfs <- mget(args, envir = environment())
  grid <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  out <- parent$FromDataAndTo("k_Desorption") |>
    dplyr::left_join(grid, by=c("Scale" = "Scale", "from.SubCompart" = "SubCompart", "Species"="Species"))|>
    dplyr::left_join(MTC_2sd, by=c("to.SubCompart" = "SubCompart"), suffix = c("_from", "_to")) |>
    
    dplyr::mutate(
      k_Desorption = dplyr::case_when(
        SpeciesName == 'Molecular' & as.character(Test) == "TRUE" & to.SubCompart == "lake" ~ NA,
        ScaleName %in% c("Regional", "Continental") & to.SubCompart == "deepocean" ~ NA,
        ScaleName %in% c("Tropic", "Moderate", "Arctic") & to.SubCompart != "deepocean" ~ NA,
        SpeciesName == 'Molecular' ~ ((MTC_2sd*MTC_2w)/(MTC_2sd + MTC_2w)/Ksdcompw ) / VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Desorption)) |>
    dplyr::select(names, k_Desorption) |>
    dplyr::arrange(dplyr::across(names))
    
  return(data.frame(out))                 
                     
  # 
  # if ((ScaleName %in% c("Regional", "Continental")) & to.SubCompartName == "deepocean") {
  #   return(NA)
  # }
  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "deepocean") {
  #   return(NA)
  # }
  # 
  # if (SpeciesName == "Molecular" & as.character(Test) == "TRUE" & to.SubCompartName == "lake"){
  #   return(NA)
  # }
  # 
  # switch (SpeciesName,
  #   "Molecular" = {
  #     # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName == "sea") {
  #     #   return(NA)
  #     # }
  #     ( (to.MTC_2sd*MTC_2w)/(to.MTC_2sd + MTC_2w)/Ksdcompw ) /
  #       VertDistance
  #   },
  #   return(NA)
  # )
  # 
}

