#' @title escape
#' @name k_Escape
#' @description calculate k for escape from air compartment to stratosphere based on t_half_Escape
#' @param t_half_Escape Half life time in air [s] 
#' @param SubCompartName considered subcompartment
#' @return k_Escape
#' @export
k_Escape <- function(t_half_Escape, parent){
  
  out <- parent$FromDataAndTo("k_Escape") |>
    dplyr::mutate(
      t_half_Escape = t_half_Escape,
      k_Escape = dplyr::case_when(
        from.SubCompart == 'cloudwater' | to.SubCompart == 'cloudwater' ~ NA,
        to.SubCompart == 'air' ~ log(2) / (t_half_Escape),
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Escape)) |>
    dplyr::select(process, from.SubCompart, to.SubCompart, Scale, Species, k_Escape) |>
    dplyr::arrange(process, from.SubCompart, to.SubCompart, Scale, Species)
  
  return(data.frame(out))
  
  # # an exclusion of cloudwater is needed as this is now also seen as an air compartment.
  # if(to.SubCompartName == "cloudwater") return (NA)
  # if(from.SubCompartName == "cloudwater") return (NA)
  # 
  # switch (to.SubCompartName,
  #   "air" =   return(log(2)/(t_half_Escape)),
  #   return(NA)
  # )

}
