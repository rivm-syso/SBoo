#'@title Run-off of particles
#'@name k_Runoff
#'@description Runoff rate corrected for soil penetration depth for molecules and particulates
#'@param RunoffFlow RunoffFlow, as computed by v_RunoffFlow [s-1]
#'@param FracROWatComp fraction of water component runoff can flow to [-]
#'@param Volume volume of compartment 
#'@param Kscompw soil water partitioning coefficient [-]
#'@param relevant_depth_s the soil depth at which a process occurs [m]
#'@param penetration_depth_s assumed penetration depth from Hollander et al. (2007), see f_CORRsoil [m]
#'@param ScaleName considered ScaleName
#'@param SubCompartName subcompartment considered
#'@param SpeciesName considered species 
#'@param Matrix type of supcompartment 
#'@param FRACrun Fraction of precipation of run-off #[-]
#'@param RAINrate Average precipitation #[m/s]
#'@param VertDistance Mixing depth soil #[m]
#'@param kinterception Fraction of particles that stays in soil - defined only for microplastics for now
#'@return k_Run-off Run-off of particles from soil to water #[s-1]
#'@export

k_Runoff <- function(RunoffFlow,VertDistance, FracROWatComp,
                     Volume, Kscompw,
                     relevant_depth_s, penetration_depth_s,
                     SpeciesName, Matrix, IntrcptFrac, parent){
  
  out <- parent$FromDataAndTo("k_Runoff")
  # Names of dataframe
  names <- dplyr::all_of(names(out))
  #convenience rename for automatic matching columns (by = ...)
  out <- dplyr::rename(out, SubCompart = from.SubCompart)
  
  # Build 
  args <- c("out", "SpeciesName", "Matrix", "RunoffFlow","VertDistance", "Volume", "Kscompw", "IntrcptFrac")
  dfs <- mget(args, envir = environment())
  out <- purrr::reduce(dfs, ~ dplyr::left_join(.x, .y)) # |> parent$states$clipStates()
  # and back
  out <- dplyr::rename(out, from.SubCompart = SubCompart) |>

    dplyr::left_join(FracROWatComp, by = c("to.SubCompart" = "SubCompart", "Species", "Scale"))

  # Output
  out <- out |>

    dplyr::mutate( 
      k_Runoff = dplyr::case_when(
        Scale %in% c("Regional", "Continental") & to.SubCompart== "sea" ~ NA,
        Scale %in% c("Tropic", "Moderate", "Arctic") & to.SubCompart!= "sea" ~ NA,
        Species == 'Unbound' ~ (RunoffFlow / Kscompw) * f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / Volume  * FracROWatComp,
        Species != 'Unbound' ~ (RunoffFlow * f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s))/ Volume * FracROWatComp * (1-IntrcptFrac), #adjusted runoff for microplastics  
        TRUE ~ NA
      )
    ) |>
    dplyr::filter(!is.na(k_Runoff)) |>
    dplyr::select(names, k_Runoff) |>
    dplyr::arrange(dplyr::across(names))
  
  return(data.frame(out))
  
  
  
  # # browser()
  # if (ScaleName %in% c("Regional", "Continental") & to.SubCompartName == "sea") {
  #   return(NA)
  # } 
  # if ((ScaleName %in% c("Tropic", "Moderate", "Arctic")) & to.SubCompartName != "sea") {
  #   return(NA)
  # } 
  # 
  # switch(SpeciesName,
  #        "Molecular" = {
  #          (RunoffFlow / Kscompw) * 
  #            f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s) / Volume  * to.FracROWatComp  #[s-1]
  #        },
  #        
  #        (RunoffFlow *
  #           f_CORRsoil(VertDistance, relevant_depth_s, penetration_depth_s))/ Volume * to.FracROWatComp * (1-IntrcptFrac) #adjusted runoff for microplastics  
  # )
}