
#' @title AddAbbreviationsSBxlsx
#' @name AddAbbreviationsSBxlsx
#' @description Returns a table, default for World$kaas with two new columns to and from with the Abbrevaition of each SimpleBox compartment
#' @param kaas table with columns toSpecies, toSubCompart and toScale, idem from...
#' @param SubcompartmentsMap names character vector with the abreviations for each SubCompart
#' @param ScalesMap names character vector with the abreviations for each Scale
#' @param SpeciesMap names character vector with the abreviations for each Species
#' @return Table with standard SimpleBox abbreviations for each compartment
#' @export
f_AddAbbreviationsSBxlsx <- function(kaas = as_tibble(World$kaas), # table with columns to... and from...
                                   SubcompartmentsMap = # all subcompartments defined for SBoo
                                     c("marinesediment" = "sd2",
                                       "freshwatersediment" = "sd1",
                                       "lakesediment" = "sd0",
                                       "agriculturalsoil" = "s2",
                                       "naturalsoil" = "s1",
                                       "othersoil" = "s3",
                                       "air" = "a",
                                       "deepocean" = "w3",
                                       "sea" = "w2",
                                       "river" = "w1",
                                       "lake" = "w0",
                                       "cloudwater" = "cw"), 
                                   ScalesMap = # all scales defined for SBoo
                                     c("Arctic" = "A",
                                       "Moderate" = "M",
                                       "Tropic" = "T",
                                       "Continental" = "C",
                                       "Regional" = "R"), 
                                   SpeciesMap = # all species defined for SBoo
                                     c("Dissolved" = "D",
                                       "Gas" = "G",
                                       "Large" = "P",
                                       "Small" = "A",
                                       "Solid" = "S",
                                       "Unbound" = "U")
) {
  
  kaas <- kaas |> mutate(
    from = paste0(
      SubcompartmentsMap[fromSubCompart],
      ScalesMap[fromScale],
      SpeciesMap[fromSpecies]
    ),
    to = paste0(
      SubcompartmentsMap[toSubCompart],
      ScalesMap[toScale],
      SpeciesMap[toSpecies]
    )
  )
  
  # kaas <-
  #   kaas |>
  #   mutate(
  #     from =
  #       ifelse((fromScale == "Tropic" | fromScale == "Arctic" | fromScale == "Moderate") &
  #                (fromSubCompart == "marinesediment" | fromSubCompart == "naturalsoil"),
  #              str_replace_all(from, c("sd2" = "sd", "s1" = "s")),
  #              from
  #       )
  #   ) |>
  #   mutate(to = ifelse((toScale == "Tropic" | toScale == "Arctic" | toScale == "Moderate") &
  #                        (toSubCompart == "marinesediment" | toSubCompart == "naturalsoil"), str_replace_all(to, c("sd2" = "sd", "s1" = "s")), to))
  # 
  return(kaas)
}