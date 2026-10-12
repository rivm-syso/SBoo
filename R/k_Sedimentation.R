#' @title k_Sedimentation
#' @name k_Sedimentation
#' @description Calculate the rate constant for sedimentation [s-1]
#' @param FRinw Fraction chemical dissolved in water [-]
#' @param VertDistance Compartment depth/height [m]
#' @param DynViscWaterStandard Dynamic viscosity of the fluid matrix
#' @param RadCP Radius of the Coarse natural particle [m]
#' @param RhoCP Density of the Coarse natural particle [m]
#' @param SettlingVelocity Settling velocity of particulate species [m.s-1]
#' @param SubCompartName Name of relevant subcompartment for which k_Sedimentation is being calculated
#' @param ScaleName Name of relevant scale for which k_Sedimentation is being calculated
#' @param SpeciesName Name of relevant species (Molecular or particulate) for which k_Sedimentation is being calculated
#' @param Test determines if SB4-Excel approach is taken or enhanced method from R version [boolean]
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_Sedimentation, the rate constant for sedimentation as first order process
#' @export
k_Sedimentation <- function(FRinw, SettlingVelocity, DynViscWaterStandard, rhoMatrix, Matrix,
                            VertDistance, RhoCP, RadCP, RadS,
                            SpeciesName, SubCompartName, ScaleName, Test, Regional_and_Continental_deepocean,parent){
  
  # Names of dataframe
  names <- dplyr::all_of(names(parent$FromDataAndTo("k_Sedimentation")))
  
  # Build grid
  args <- c("SpeciesName", "SubCompartName", "ScaleName", "Matrix", "FRinw", "SettlingVelocity", "rhoMatrix", "VertDistance", "RhoCP", "RadCP")
  dfs <- mget(args, envir = environment())
  grid <- purrr::reduce(dfs, ~ dplyr::full_join(.x, .y, by= intersect(names(.x), names(.y)))) |> parent$states$clipStates()
  
  # Output
  out <- parent$FromDataAndTo("k_Sedimentation") |>
    # From
    dplyr::left_join(grid, by=c("from.SubCompart" = "SubCompart", "Scale"="Scale", "Species"="Species")) |>
    
    dplyr::mutate(
      k_Sedimentation = dplyr::case_when(
        (Scale %in% c("Tropic", "Moderate", "Arctic")) & from.SubCompart == "sea" & to.SubCompart== "marinesediment" ~ 
          NA,
        (Scale %in% c("Regional", "Continental")) &
          ((from.SubCompart == "deepocean" & to.SubCompart == "marinesediment") | (from.SubCompart == "sea" & to.SubCompart == "deepocean")) &
          (isFALSE(Regional_and_Continental_deepocean) | is.na(Regional_and_Continental_deepocean) | Regional_and_Continental_deepocean == "FALSE") ~
          NA,
        (Scale %in% c("Regional", "Continental")) & from.SubCompart == "sea" & to.SubCompart == "marinesediment" &
          (Regional_and_Continental_deepocean == "TRUE" | isTRUE(Regional_and_Continental_deepocean)) ~
          NA,
        (Scale %in% c("Tropic", "Moderate", "Arctic")) & from.SubCompart %in% c("lake","river") ~ NA,
        SpeciesName == 'Molecular' & to.SubCompart == 'deepocean' ~ NA,
        SpeciesName == 'Molecular' & as.character(Test) == "TRUE" & to.SubCompart == 'lakesediment' ~ NA,
        SpeciesName == 'Molecular' & as.character(Test) == "TRUE" & to.SubCompart != 'lakesediment' ~
          (2.5/(24*3600))*(1 - FRinw) / VertDistance,
        SpeciesName == 'Molecular' & as.character(Test) != "TRUE" ~
          (f_SetVelWater(Shortest_side=RadCP*2, rho_species=RhoCP, rhoMatrix=rhoMatrix, 
                        DynViscWaterStandard=DynViscWaterStandard,DynViscAirStandard=NA,
                        Matrix=Matrix,SubCompartName=from.SubCompart, Shape=NA,Longest_side=NA, 
                        Intermediate_side=NA, DragMethod="Original")) * (1 - FRinw) / VertDistance,
        SpeciesName != 'Molecular' & (SettlingVelocity <= 0 | is.na(SettlingVelocity)) ~ 0,
        SpeciesName != 'Molecular' ~ SettlingVelocity/VertDistance,
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Sedimentation)) |>
    dplyr::select(names, k_Sedimentation) |>
    dplyr::arrange(dplyr::across(names))
  
  return(data.frame(out))
  
  # 
  # # Sedimentation at global scales goes from sea to deeopocean to marinesediment
  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName == "sea" & to.SubCompartName == "marinesediment") {
  #   return(NA)
  # }
  # 
  # # If Regional_and_Continental_deepocean is TRUE, the sedimentation rate from sea to deepocean and deepocean to marinesediment
  # # at Regional and Continental scale should be NA, because the sedimentation rate goes directly from sea to marinesediment
  # # at these scales.  
  # if ((ScaleName %in% c("Regional", "Continental")) &&
  #     ((SubCompartName == "deepocean" && to.SubCompartName == "marinesediment") || (SubCompartName == "sea" && to.SubCompartName == "deepocean")) &&
  #     (isFALSE(Regional_and_Continental_deepocean) || is.na(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE")) {
  #   return(NA)
  # }
  # 
  # # If Regional_and_Continental_deepocean is TRUE, remove the sedimentation rate from sea to marinesediment at Regional and Continental scale
  # # because deepocean is added between sea and deeopocean.
  # if ((ScaleName %in% c("Regional", "Continental")) && SubCompartName == "sea" && to.SubCompartName == "marinesediment" &&
  #     (Regional_and_Continental_deepocean == "TRUE" || isTRUE(Regional_and_Continental_deepocean))) {
  #   return(NA)
  # }
  # 
  # # No lake and river present at global scales
  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & SubCompartName %in% c("lake","river")) {
  #   return(NA)
  # }
  # 
  # switch(SpeciesName,
  #        "Molecular" = {
  #          if (to.SubCompartName == "deepocean") {
  #            return(NA)
  #          } 
  #          if (as.character(Test) == "TRUE") {
  #            if (to.SubCompartName == "lakesediment") {
  #              return(NA)
  #            } else {
  #              SetlingVelocityCP <- 2.5/(24*3600)
  #              return(SetlingVelocityCP*(1 - FRinw) / VertDistance)
  #            }
  #          }
  #          SetlingVelocityCP <- 
  #            f_SetVelWater(Shortest_side=from.RadCP*2, 
  #                          rho_species=from.RhoCP, 
  #                          rhoMatrix=rhoMatrix, 
  #                          DynViscWaterStandard=DynViscWaterStandard,
  #                          DynViscAirStandard=NA,
  #                          Matrix=Matrix,SubCompartName=SubCompartName, 
  #                          Shape=NA,
  #                          Longest_side=NA, Intermediate_side=NA,
  #                          DragMethod="Original")
  #          
  #          return(SetlingVelocityCP*(1 - FRinw) / VertDistance)
  #        },
  #        { 
  #          if (SettlingVelocity <= 0 | is.na(SettlingVelocity)) { # The is.na part was added because for some particles when MinSettVel is na the setvel is not recalculated for some reason which results in an NA value and thus an error in the if statement. 
  #            return(0)
  #          }
  #          return(SettlingVelocity/VertDistance)
  #        }
  #        
  # )
  
}