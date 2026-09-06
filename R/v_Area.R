#' @title Area
#' @name Area
#' @param AREAland the land area in the considered compartment [m2]
#' @param AREAsea  the sea area in the considered compartment [m2]
#' @param landFRAC fraction of land in the considered compartment [-]
#' @param SubCompartName the subcompartment which is considered
#' @param ScaleName the scale which is considered
#' @return Area (but not for sediment) based on data for the SubCompartment / Scale
#' @export
Area <- function (AreaLand,
                  AreaSea,
                  landFRAC,
                  SubCompartName,
                  ScaleName,
                  Regional_and_Continental_deepocean,
                  parent,
                  SpeciesName) {
  
  out <- ScaleName |>
    expand_grid(SubCompartName, SpeciesName) |>
    full_join(AreaSea, "Scale") |>
    full_join(AreaLand, "Scale") |>
    full_join(landFRAC, c("Scale", "SubCompart")) |>
    parent$states$clipStates() |>
    group_by(ScaleName) |>
    mutate(
      landFRAC_used = dplyr::case_when(
        SubCompartName == "lakesediment" &
          ScaleName %in% c("Regional", "Continental") ~
          landFRAC[SubCompart == "lake"][1],
        
        SubCompartName == "freshwatersediment" &
          ScaleName %in% c("Regional", "Continental") ~
          landFRAC[SubCompart == "river"][1],
        
        SubCompartName %in% c("agriculturalsoil", "naturalsoil",
                              "othersoil", "lake", "river") ~ landFRAC,
        
      ),
      Area = dplyr::case_when(
        SubCompartName %in% c("air", "cloudwater") ~ AreaLand + AreaSea,

        SubCompartName %in% c("sea", "marinesediment") ~ AreaSea,

        SubCompartName %in% c("deepocean") & ScaleName %in% c("Arctic", "Moderate", "Tropic") ~ AreaSea,

        SubCompartName %in% c("deepocean") &
          ScaleName %in% c("Regional", "Continental") &
          (!is.na(Regional_and_Continental_deepocean) &&
             (isTRUE(Regional_and_Continental_deepocean) ||
                Regional_and_Continental_deepocean == "TRUE")) ~ AreaSea,

        SubCompartName %in% c(
          "lakesediment", "freshwatersediment",
          "agriculturalsoil", "othersoil", "naturalsoil", "lake", "river"
        ) &
          ScaleName %in% c("Regional", "Continental") ~
          landFRAC_used * AreaLand,

        !(ScaleName %in% c("Regional", "Continental")) & SubCompartName %in% c('othersoil', 'naturalsoil') ~ AreaLand, 
        TRUE ~ NA_real_
      )
    ) |>
    ungroup() |>
    select(Scale, SubCompart, Area)  
  
  # Add area for cloudwater
  extra <- out |>
    filter(SubCompart == 'air') |>
    mutate(
      SubCompart = 'cloudwater'
    )
  out <- out |>
    bind_rows(extra) |>
    filter(!is.na(Area)) |>
    arrange(Scale, SubCompart, Area)
  
  return(data.frame(out))

  # # easiest
  # if (SubCompartName %in% c("air", "cloudwater")) {
  #   return(AreaLand + AreaSea)
  # }
  # 
  # if (SubCompartName %in% c("sea", "marinesediment")) {
  #   return(AreaSea)
  # }
  # if (SubCompartName == "deepocean" &
  #     ScaleName %in% c("Arctic", "Moderate", "Tropic")) {
  #   return(AreaSea)
  # }
  # if (SubCompartName == "deepocean" &
  #     ScaleName %in% c("Regional", "Continental") &
  #     (!is.na(Regional_and_Continental_deepocean) && (isTRUE(Regional_and_Continental_deepocean) || Regional_and_Continental_deepocean == "TRUE"))) {
  #   return(AreaSea)
  # }
  # if (SubCompartName == "lakesediment" & ScaleName %in% c("Regional", "Continental")){
  #   return(all.landFRAC$landFRAC[all.landFRAC$SubCompart == "lake" & all.landFRAC$Scale == ScaleName] *AreaLand)
  # }
  # 
  # if (SubCompartName == "freshwatersediment" & ScaleName %in% c("Regional", "Continental")){
  #   return(all.landFRAC$landFRAC[all.landFRAC$SubCompart == "river" & all.landFRAC$Scale == ScaleName] *AreaLand)
  # }
  # 
  # 
  # # on land, lake, freshwater and soils;
  # if (ScaleName %in% c("Regional", "Continental")) {
  #   return(landFRAC * AreaLand)
  # }
  # 
  # # on land, on global scales: only naturalsoil 
  # #! Bij deze IF zijn alleen globale scales nog over, de parent clip states clipped natural soil uit deze schalen, 
  # #hier is enkel othersoil beschikbaar.
  # if (SubCompartName == "naturalsoil") {
  #   return(AreaLand)
  # }
  # 
  # #all other cases
  # return(NA)
}
