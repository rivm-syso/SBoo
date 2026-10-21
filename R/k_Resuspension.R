#' @title Resuspension rate constant of substances in sediment
#' @name k_Resuspension
#' @description Calculation of the resuspension rate. The top layer of sediment is assumed to be well mixed and thus continiously refreshed,
#' for additional details see Schoorl et al. (2015)
#' @param DynViscWaterStandard Dynamic viscosity of water [kg m-1 s-1]
#' @param rhoMatrix density of the Matrix [kg m-3]
#' @param NETsedrate net sedimentation rate, data input [s] (Schoorl et al., 2015)
#' @param VertDistance mixed depth sediment compartment #[m]
#' @param RhoCP Mineral density of sediment and soil #[kg/m3]
#' @param FRACs Volume fraction solids in sediment #[-]
#' @param RadCP radius of coarse particulate particles [m]
#' @param SUSP mass concentration of suspended matter [kg m-3]
#' @param SpeciesName species considered
#' @param ScaleName scale considered
#' @param SubCompartName subcompartment considered
#' @param Test determines if SB4-Excel approach is taken or enhanced method from R version [boolean]
#' @param SettlingVelocitySPM settling velocity of suspended matter particles
#' @param Regional_and_Continental_deepocean If this variable is TRUE, Regional and Continental deepocean compartments are removed
#' @return k_Resuspension Resuspension flow from sediment #[s-1]
#' @export

k_Resuspension <- function(VertDistance, # SettlVelocitywater
                           DynViscWaterStandard,
                           rhoMatrix,
                           NETsedrate,
                           RadCP, RhoCP, FRACs, SUSP, 
                           Matrix,
                           SpeciesName, ScaleName, SubCompartName, Test, Regional_and_Continental_deepocean,
                           SettlingVelocity, parent) {
  # Making to dataframe
  to_water <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    dplyr::full_join(Matrix, by ="SubCompart") |>
    parent$states$clipStates() |>
    dplyr::rename(to.Matrix = Matrix) |>
    dplyr::filter(to.Matrix %in% c("water") & SubCompart != 'cloudwater') |>
    dplyr::rename(to.SubCompart = SubCompart) |>
    dplyr::mutate(
      SubCompart = dplyr::case_when(
        to.SubCompart == 'deepocean' ~ 'marinesediment',
        to.SubCompart == 'sea' ~ 'marinesediment',
        to.SubCompart == 'river' ~ 'freshwatersediment',
        to.SubCompart == 'lake' ~ 'lakesediment',
        TRUE ~ NA
      )
    ) |>
    # dplyr::left_join(RadCP, by = c("to.SubCompart" = "SubCompart")) |>
    # dplyr::left_join(RhoCP, by = c("to.SubCompart" = "SubCompart")) |>
    # dplyr::left_join(rhoMatrix, by = c("to.Matrix" = "Matrix")) |>
    dplyr::left_join(SUSP, by = c("to.SubCompart" = "SubCompart"))|>
    dplyr::left_join(NETsedrate, by = c("to.SubCompart" = "SubCompart", "Scale" = "Scale")) |>
    dplyr::left_join(SettlingVelocity, by = c("Scale"= "Scale", "to.SubCompart" = "SubCompart", "Species" = "Species")) |>
    dplyr::select(-SubCompartName, -SpeciesName, -ScaleName)

  # Making from dataframe
  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    dplyr::full_join(VertDistance) |>
    dplyr::full_join(Matrix) |>
    dplyr::full_join(FRACs) |>
    dplyr::full_join(RhoCP) |>
    parent$states$clipStates() |>
    dplyr::filter(Matrix == 'sediment') |>
    dplyr::left_join(to_water, by=c("SubCompart" = "SubCompart", "Scale" = "Scale", "Species" = "Species")) |>
    dplyr::mutate(
      SettlingVelocity = dplyr::case_when(
        SpeciesName == "Molecular" & as.character(Test) == "TRUE" ~ 2.5 / (24 * 3600),
        TRUE ~ SettlingVelocity
      ),
      

      # Gross sedimentation rate from water [m/s]
      GROSSEDrate = SettlingVelocity * SUSP / (FRACs * RhoCP),
      
      # Resuspension flow from sediment [m/s]; can't be < 0
      RESUSflow = pmax(0, GROSSEDrate - NETsedrate),

      # Resuspension k to water [s-1]
      k_Resuspension = dplyr::case_when(
        # If Regional_and_Continental_deepocean is FALSE, no resuspension from marinesediment to deepocean
        ScaleName %in% c("Regional", "Continental") & to.SubCompart == 'deepocean' &
          (isFALSE(Regional_and_Continental_deepocean) || is.na(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE") ~
        NA,

        #If Regional_and_Continental_deepocean is TRUE, no resuspension from marinesediment to sea
        ScaleName %in% c("Regional", "Continental") & to.SubCompart == 'sea' &
          (!is.na(Regional_and_Continental_deepocean) && (isTRUE(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "TRUE")) ~
        NA,
        TRUE ~ RESUSflow / VertDistance
      )
    ) |>
    dplyr::rename(from.SubCompart = SubCompart) |>
    dplyr::filter(!is.na(k_Resuspension)) |>
    dplyr::arrange(Scale, from.SubCompart, to.SubCompart, Species) |>
    dplyr::select(Scale, from.SubCompart, to.SubCompart, Species, k_Resuspension)

  return(data.frame(out))
  
  # 
  # 
  # 
  # # If Regional_and_Continental_deepocean is FALSE, no resuspension from marinesediment to deepocean
  # if ((ScaleName %in% c("Regional", "Continental")) &&
  #          (to.SubCompartName == "deepocean") &&
  #          (isFALSE(Regional_and_Continental_deepocean) || is.na(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "FALSE")) {
  #   return(NA)
  # }
  # 
  # #If Regional_and_Continental_deepocean is TRUE, no resuspension from marinesediment to sea
  # if ((ScaleName %in% c("Regional", "Continental")) &&
  #     (to.SubCompartName == "sea") &&
  #     (!is.na(Regional_and_Continental_deepocean) && (isTRUE(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "TRUE"))
  # ) {
  #   return(NA)
  # } 
  # 
  # else if (SpeciesName == "Molecular") {
  #   if (as.character(Test) == "TRUE") {
  #     SettlingVelocitySPM <- 2.5 / (24 * 3600)
  #   } else {
  #     # ScaleName
  #     SettlingVelocitySPM <- 
  #       f_SetVelWater(Shortest_side=to.RadCP*2, 
  #                     rho_species=to.RhoCP, 
  #                     rhoMatrix=to.rhoMatrix, 
  #                     DynViscWaterStandard=DynViscWaterStandard,
  #                     DynViscAirStandard=NA,
  #                     Matrix=to.Matrix,SubCompartName=to.SubCompartName, 
  #                     Shape=NA,
  #                     Longest_side=NA, Intermediate_side=NA,
  #                     DragMethod="Original")
  #     
  #   }
  # } else {
  #   # ScaleName
  #   SettlingVelocitySPM <-  
  #     f_SetVelWater(Shortest_side=to.RadCP*2, 
  #                   rho_species=to.RhoCP, 
  #                   rhoMatrix=to.rhoMatrix, 
  #                   DynViscWaterStandard=DynViscWaterStandard,
  #                   DynViscAirStandard=NA,
  #                   Matrix=to.Matrix,SubCompartName=to.SubCompartName, 
  #                   Shape=NA,
  #                   Longest_side=NA, Intermediate_side=NA,
  #                   DragMethod="Original")
  # }
  # 
  # # Gross sedimentation rate from water [m/s]
  # GROSSEDrate <- SettlingVelocitySPM * to.SUSP / (FRACs * from.RhoCP) # [m.s-1] possibly < NETsedrate
  # 
  # # Resuspension flow from sediment [m/s]; can't be < 0
  # RESUSflow <- max(0, GROSSEDrate - to.NETsedrate) # for particulates this NETsedrate is not optimal!
  # 
  # # Resuspension k to water [s-1]
  # return(RESUSflow / VertDistance)
}
