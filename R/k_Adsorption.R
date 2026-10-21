#'@title Adsorption of chemical to across an interface
#'@name k_Adsorption
#'@description The adsorption rate constant (k_Adsorption) is calculated based on the principles of gas absorption as described by Hollander
#'et al. (2004) 10.1080/10629360412331297470 & Hollander et al. (2007) https://doi.org/10.1016/j.envpol.2006.09.018
#'@param FRinw fraction of species in water,see FRinw()
#'@param MTC_2sd mass transfer coefficient to sediment, see MTC_2sd()
#'@param MTC_2w mass transfer coefficient to water, see SBvariable
#'@param MTC_2a mass transfer coefficient to air, see SBvariable
#'@param MTC_2s mass transfer coefficient to soil, see SBvariable
#'@param FRorig fraction of species in original form, see SBvariable
#'@param Kacompw partitioning coefficient for air, see SBvariable
#'@param Kscompw partitioning coefficient for soil, see SB variable 
#'@param Matrix Receiving compartment type 
#'@param to.Area Area of the receiving compartent
#'@param VertDistance Vertical distance of the compartment 
#'@param AreaLand calculated land area of the scale considered, see SBvariable
#'@param AreaSea calculated sea area of the scale considered, see SBvariable 
#'@param FRorig_spw fraction of original species in soil pore water 
#'@param SubCompartName name of the subcompartname that is from/to 
#'@param ScaleName name of scale
#'@param Test Test = TRUE or FALSE, depending on if SBExcel version should be copied
#'@param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#'@returns The adsorption rate constant relevant for the receiving compartments soil, water or sediment [s-1]
#'@export
k_Adsorption <- function (FRingas, FRinw, MTC_2sd, FRorig_spw,
                          MTC_2w, MTC_2a, MTC_2s, FRorig, Kacompw, 
                          Kscompw, Matrix, VertDistance, 
                          AreaLand, AreaSea, Area,
                          SubCompartName, ScaleName, Test, Regional_and_Continental_deepocean,
                          parent, SpeciesName) {
  
  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Adsorption")))
  
  # copies of whole dataframes for the test argument
  FRorig_og = FRorig
  FRorig_spw_og = FRorig_spw
  
  # Build grid from
  args <- c("SubCompartName", "ScaleName", "SpeciesName", "Matrix","FRingas", "FRinw", "MTC_2sd","MTC_2w", "MTC_2s", "Kacompw", "VertDistance", "AreaLand", "AreaSea")
  dfs <- mget(args, envir = environment())
  grid_from <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  # Build grid to
  args <- c("SubCompartName", "ScaleName", "SpeciesName", "Matrix","FRorig", "Kscompw", "FRorig_spw","MTC_2a", "Area", "MTC_2w")
  dfs <- mget(args, envir = environment())
  grid_to <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  out <- parent$FromDataAndTo("k_Adsorption") |>
    # From
    dplyr::left_join(grid_from, by=c("Scale" = "Scale", "from.SubCompart" = "SubCompart", "Species"= "Species")) |>
    
    # To met suffix 
    dplyr::left_join(grid_to, by=c("Scale" = "Scale", "to.SubCompart" = "SubCompart", "Species"= "Species"), suffix = c("_from", "_to")) |>
    
    dplyr::mutate(
      FRorig = dplyr::case_when(
        Matrix_to == 'water' & as.character(Test) == 'TRUE' ~ dplyr::filter(FRorig_og, SubCompart =='river') |> dplyr::pull(FRorig),
        TRUE ~ FRorig
      ),
      FRorig_spw = dplyr::case_when(
        Matrix_to == 'soil' & as.character(Test) == 'TRUE' ~ dplyr::filter(FRorig_spw_og, SubCompart =='naturalsoil') |> dplyr::pull(FRorig_spw),
        TRUE ~ FRorig_spw
      ),
      GASABS = dplyr::case_when(
        Matrix_to == 'water' ~ FRingas*(MTC_2w_from*MTC_2a/(MTC_2w_from*(Kacompw*FRorig)+MTC_2a)),
        Matrix_to == 'soil' ~ FRingas*(MTC_2s*MTC_2a)/(MTC_2s*(Kacompw*FRorig_spw)/Kscompw+MTC_2a),
        TRUE ~ NA
      ),
      AreaFrac = dplyr::case_when(
        Matrix_to == 'water' ~ Area/(AreaLand+AreaSea),
        Matrix_to == 'soil' ~ Area/(AreaLand+AreaSea),
        TRUE ~ NA
      ),
      ADSORB = dplyr::case_when(
        Matrix_to == 'sediment' ~ (MTC_2sd*MTC_2w_to)/(MTC_2sd+MTC_2w_to)*FRinw,
        TRUE ~ NA
      ),
      k_Adsorption = dplyr::case_when(
        Matrix_to == 'water' ~ GASABS/VertDistance*AreaFrac,
        Matrix_to == 'soil' ~ GASABS/VertDistance*AreaFrac,
        Matrix_to == 'sediment' & as.character(Test) == "TRUE" & to.SubCompart == "lakesediment" ~ NA,
        Matrix_to == 'sediment' ~ ADSORB/VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Adsorption)) |>
    dplyr::select(names, k_Adsorption) |>
    dplyr::arrange(dplyr::across(names))

  return(data.frame(out))

  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & from.SubCompartName == "sea") {
  #   return(NA)
  # }
  # if ((isFALSE(Regional_and_Continental_deepocean) || is.na(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE") &&
  #   (ScaleName %in% c("Regional", "Continental") & from.SubCompartName == "deepocean")) {
  #   return(NA)
  # }
  # switch(to.Matrix,
  #        
  #        "water" = { # air to water
  #          if (ScaleName %in% c("Regional", "Continental")){
  #            if (as.character(Test) == "TRUE"){
  #              to.FRorig <-  all.FRorig |>
  #                filter(SubCompart == "river") 
  #              to.FRorig <- to.FRorig$FRorig
  #            } 
  #          }
  #          GASABS = FRingas*(from.MTC_2w*to.MTC_2a/(from.MTC_2w*(Kacompw*to.FRorig)+to.MTC_2a))
  #          AreaFrac = to.Area/(AreaLand+AreaSea)
  #          return(GASABS/VertDistance*AreaFrac) },
  #        "soil" = { # air to soil
  #          if (as.character(Test) == "TRUE"){
  #            to.FRorig_spw <- all.FRorig_spw |>
  #              filter(SubCompart == "naturalsoil")
  #            to.FRorig_spw <- to.FRorig_spw$FRorig_spw
  #          } 
  #          GASABS = FRingas*(from.MTC_2s*to.MTC_2a)/(from.MTC_2s*(Kacompw*to.FRorig_spw)/to.Kscompw+to.MTC_2a)
  #          AreaFrac = to.Area/(AreaLand+AreaSea)
  #          return(GASABS/VertDistance*AreaFrac) },
  #        "sediment" = { # water to sediment
  #          ADSORB = (from.MTC_2sd*to.MTC_2w)/(from.MTC_2sd+to.MTC_2w)*FRinw
  #          if (as.character(Test) == "TRUE" && to.SubCompartName == "lakesediment"){
  #            return(NA)
  #          } else {
  #             return(ADSORB/VertDistance) 
  #          }
  #         }, 
  #        return(NA)
  # )
  
}
