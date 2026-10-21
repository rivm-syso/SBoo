#' @title from flow to Advection "process"
#' @name k2_Advection
#' @description Calculation of k, given a Flow-dataframe
#' @param flow advection rate [m/s]
#' @param Volume volume of compartment [m3]
#' @return Rate constant for 1st order process associated with fluxes
#' @importFrom rlang :=
#' @importFrom rlang !!
#' @importFrom rlang sym
#' @export
#' 

k2_Advection <- function(flow, Volume) {
  
  # dimension that were not the changing dimension 
  dimNames <- colnames(flow)[colnames(flow) %in% The3D]
  # the dimension that was
  dimChanged <- The3D[! The3D %in% dimNames]
  # easier to match for the dplyr::leftjoin
  toChanged <- paste0("to.", dimChanged)
  
  out <- flow |>
    dplyr::mutate(!!dimChanged := .data[[toChanged]]) |>  # copy values from to.age into age
    dplyr::left_join(Volume) |>            # 
    dplyr::mutate(k = flow / Volume) |>
    dplyr::select(-dplyr::all_of(dimChanged), -Volume)  

  # expand to from for all dimensions
  toNames <- paste0("to.", dimNames)
  fromNames <- paste0("from.", dimNames)
  for (i_nm in seq_along(dimNames)) {
    out <- out |>
      dplyr::mutate(
        !!toNames[i_nm] := .data[[dimNames[i_nm]]] ) |>
      dplyr::rename(!!fromNames[i_nm] := !!dimNames[i_nm])
  }

  return(out) 
}