#'@title Burial by sediment
#'@name k_Burial
#'@description Burial based on net sedimentation rates, as described in Schoorl et al. (2015)
#'@param VertDistance mixed depth water sediment compartment [m]
#'@param NETsedrate Net sediment accumulation rate (from the surface water above). Values are constants as reported in Schoorl et al. (2015) [m/s]
#'@param ScaleName Scale name of the considered process
#'@param SubCompartName Subcompartment name for which the computation is done. 
#'@return k_Burial Burial from sediment [s-1]
#'@export

k_Burial <- function(VertDistance, NETsedrate, ScaleName, SubCompartName, parent, SpeciesName){
  
  # Get NETsedrates of the water column above the sediment, ie map&link the correct subcompartments
  NETsedrate_sediments <- NETsedrate |>
    dplyr::filter(SubCompart %in% c("lake", 'deepocean', 'sea', 'river')) |>
    dplyr::mutate(
      from.SubCompart = dplyr::case_when(
        SubCompart == 'lake' ~ 'lakesediment',
        SubCompart == 'deepocean' & Scale %in% c("Tropic", "Moderate", "Arctic") ~ 'marinesediment',
        SubCompart == 'river' ~ 'freshwatersediment',
        SubCompart == 'sea' ~ 'marinesediment',
        TRUE ~ NA
      )
    ) |>
    dplyr::rename(to.SubCompart = SubCompart) |>
    dplyr::filter(!is.na(from.SubCompart)) |> # NETsedrate has deepocean in regional whereas clip states doesnt 
    dplyr::select(NETsedrate, from.SubCompart, Scale, to.SubCompart)

  out <- ScaleName |>
    tidyr::expand_grid(SubCompartName, SpeciesName) |>
    parent$states$clipStates() |>
    dplyr::left_join(VertDistance, by=c("Scale", "SubCompart")) |>
    dplyr::filter(SubCompart %in% c("lakesediment", "marinesediment", "freshwatersediment")) |>
    dplyr::left_join(NETsedrate_sediments, by=c("Scale" = "Scale", "SubCompart" = "from.SubCompart")) |>
    dplyr::mutate(
      k_Burial = NETsedrate / VertDistance
    ) |>
    dplyr::filter(!is.na(k_Burial)) |>
    dplyr::rename(from.SubCompart = SubCompart) |>
    dplyr::arrange(Scale, from.SubCompart, to.SubCompart, Species) |>
    dplyr::select(Scale, from.SubCompart, to.SubCompart, Species, k_Burial)
    
  return(data.frame(out))

  
  # # NETsedrate assumed identical to NETsedrate of the water column above the sediment
  # waterabove <- switch (SubCompartName,
  #     "lakesediment" = "lake",
  #     "marinesediment" = {switch(ScaleName,
  #                                "Tropic" = "deepocean",
  #                                "Moderate" = "deepocean",
  #                                "Arctic" = "deepocean",
  #                                "sea")}
  #       ,
  #     "freshwatersediment" = "river",
  #     NA
  # )
  # if (is.na(waterabove)) return (NA)
  # rightRow <- which(all.NETsedrate$SubCompart == waterabove & all.NETsedrate$Scale == ScaleName)
  # if (length(rightRow) != 1) return(NA)
  # waterNETsedrate <- all.NETsedrate$NETsedrate[rightRow]
  # waterNETsedrate / VertDistance
}
