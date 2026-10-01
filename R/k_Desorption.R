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
k_Desorption <- function (Ksdcompw, MTC_2w, to.MTC_2sd, VertDistance,
                          SpeciesName, to.SubCompartName, SubCompartName, ScaleName, Test, parent) {
  
  
  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    parent$states$clipStates() |>
    dplyr::filter(SubCompartName %in% c("lakesediment", "freshwatersediment", "marinesediment")) |>
    dplyr::left_join(VertDistance, by=c("Scale", "SubCompart"))|>
    dplyr::left_join(Ksdcompw, by=c("Scale", "SubCompart")) |>
    dplyr::left_join(MTC_2w, by=c("Scale", "SubCompart")) |>
    dplyr::mutate(
      to.SubCompartName = dplyr::case_when(
        SubCompart == 'freshwatersediment' ~ 'river',
        SubCompart == 'lakesediment' ~ 'lake',
        SubCompart == 'marinesediment' & !ScaleName %in% c("Regional", "Continental")  ~ 'sea',
        SubCompart == 'marinesediment' ~ 'deepocean',
        TRUE ~NA
      ),
    ) |>
    dplyr::left_join(to.MTC_2sd, by =c("to.SubCompartName" = "SubCompart")) |>
    dplyr::mutate(
      k_Desorption = dplyr::case_when(
        SpeciesName == 'Molecular' & as.character(Test) == "TRUE" & to.SubCompartName == "lake" ~ NA,
        SpeciesName == 'Molecular' ~ ((MTC_2sd*MTC_2w)/(MTC_2sd + MTC_2w)/Ksdcompw ) / VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::select(Scale,SubCompart,Species,k_Desorption) |>
    dplyr::arrange(Scale,SubCompart,Species)
                     
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

