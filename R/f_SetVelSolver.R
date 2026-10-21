#' @title Settling Velocity Solver based on RSS
#' @name f_SetVelSolver
#' @description Calculates the settling velocity by minimizing  the residual sum of squares (RSS)
#' @param CD Drag Coefficient of a particle [-]
#' @param DragMethod Method of calculating the Drag Coefficient
#' @param Psi Shape factor, circularity/sphericity [-]
#' @param Re Reynolds number, as returned by the solver [-]
#' @param CSF Corey Shape Factor [-]
#' @param d_eq Equivalent spherical diameter of the particle [-]
#' @param DynViscFluidStandard Dynamic viscosity of liquid  (water or air) compartment [kg.m-1.s-1]
#' @param rhoParticle Density of nanoparticle [kg.m-3]
#' @param rhoFluid Density of water or air [kg m-3]
#' @param rad_species radius of the species [m]
#' @param GN gravitational force constant [m2 s-1]
#' @param Cunningham Cunningham coefficient, see f_Cunningham [-]
#' @return settling velocity [m/s] 
#' @export
#' 
f_SetVelSolver <- function(d_eq, Psi, DynViscFluidStandard, rhoParticle, rhoFluid, DragMethod, CSF, Matrix, rad_species, kS, kN) {
  
  GN <- constants::syms$gn
  n  <- length(d_eq)
  tol = 1e-10
  
  # basischecks
  stopifnot(
    length(Psi) == n,
    length(DynViscFluidStandard) == n,
    length(rhoParticle) == n,
    length(rhoFluid) == n,
    length(DragMethod) == n,
    length(CSF) == n,
    length(Matrix) == n,
    length(rad_species) == n,
    length(kS) == n,
    length(kN) == n
  )
  # d_eq = out$d_eq
  # Psi = out$Psi
  # DynViscFluidStandard = out$DynViscWaterStandard
  # rhoParticle = out$rho_species
  # rhoFluid = out$rhoMatrix
  # DragMethod = out$DragMethod
  # CSF= out$CSF
  # Matrix = out$Matrix
  # rad_species = out$rad_species
  # kS= out$kS
  # kN = out$kN
  
  # beginwaarden
  v_s <- rep(1e-4, n)
  
  # Replace indices of values where DragMethod is incompatible
  DragMethod[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  if (all(is.na(DragMethod))) {
    # Return early
    return(rep(NA, n))
  }
  d_eq[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  Psi[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  DynViscFluidStandard[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  rhoParticle[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  rhoFluid[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  CSF[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  Matrix[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  rad_species[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  kS[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  kN[!DragMethod %in% c("Dioguardi", "Default", "Swamee", "Stokes", "Bagheri")] <- NA
  
  Value_DragMethod = unique(DragMethod)[1]
  
  for (iter in seq_len(100)) {
    Re <- d_eq * v_s * rhoFluid / DynViscFluidStandard
    CD <- f_DragCoefficient(Value_DragMethod, Re, Psi, CSF, kS, kN)

    Cunningham <- ifelse(Matrix == "air", f_Cunningham(rad_species), 1)

    base_term <- sqrt(4 / 3 * d_eq / CD *
                        ((rhoParticle - rhoFluid) / rhoFluid) * GN)

    v_s_new <- base_term * Cunningham

    # simpele validatie
    v_s_new[!is.finite(v_s_new)] <- NA_real_

    diff <- abs(v_s_new - v_s)
    max_diff <- suppressWarnings(max(diff, na.rm = TRUE))
    if (!is.infinite(max_diff) && max_diff < tol) {
      v_s <- v_s_new
      break
    }

    v_s <- v_s_new
  }

  v_s
  
 
  
  # Define the RSS function to be minimized
  # GN <- constants::syms$gn
  # 
  # 
  # out <- data.frame(d_eq=d_eq, Psi=Psi, DynViscFluidStandard=DynViscFluidStandard, rhoParticle=rhoParticle, rhoFluid=rhoFluid, DragMethod=DragMethod, CSF=CSF, Matrix=Matrix, rad_species=rad_species, kS=kS, kN=kN) |>
  #   dplyr::rowwise() |>
  #   dplyr::mutate(
  #     set_vel = list(
  #       if (Matrix == "water") {
  #         RSS_function <- function(v_s) {
  #           Re <- d_eq * v_s * rhoFluid / DynViscFluidStandard
  #           CD <- f_DragCoefficient(DragMethod, Re, Psi, CSF, kS, kN)
  #           v_s_new <- sqrt(4 / 3 * d_eq / CD *
  #                             ((rhoParticle - rhoFluid) / rhoFluid) * GN)
  #           (v_s - v_s_new)^2
  #         }
  #         optimize(RSS_function, interval = c(0, 1), tol = 1e-10)$minimum
  #       } else if (Matrix == "air") {
  #         RSS_function <- function(v_s) {
  #           Re <- d_eq * v_s * rhoFluid / DynViscFluidStandard
  #           Cunningham <- f_Cunningham(rad_species)
  #           CD <- f_DragCoefficient(DragMethod, Re, Psi, CSF, kS, kN)
  #           v_s_new <- sqrt(4 / 3 * d_eq / CD *
  #                             ((rhoParticle - rhoFluid) / rhoFluid) * GN) *
  #             Cunningham
  #           (v_s - v_s_new)^2
  #         }
  #         optimize(RSS_function, interval = c(0, 10), tol = 1e-10)$minimum
  #       } else {
  #         NA_real_
  #       }
  #     )
  #   ) |>
  #   dplyr::ungroup()
  # 
  # return(out$set_vel)
  
  
  
  # Single value
  # switch(Matrix,
  #         "water" = {
  # 
  #           RSS_function <- function(v_s) {
  #             Re <- d_eq * v_s * rhoFluid / DynViscFluidStandard
  #             CD <- f_DragCoefficient(DragMethod, Re, Psi, CSF, kS, kN)
  #             v_s_new <- sqrt(4 / 3 * d_eq / CD * ((rhoParticle - rhoFluid) / rhoFluid) * GN)
  #             RSS <- (v_s - v_s_new) ^ 2
  #             return(RSS)}
  #           result <- optimize(RSS_function, interval = c(0, 1), tol = 1e-10)
  # 
  #           return(result$minimum)
  #           },
  #         "air" = {
  #           RSS_function <- function(v_s) {
  #             Re <- d_eq * v_s * rhoFluid / DynViscFluidStandard
  #             Cunningham <- f_Cunningham(rad_species)
  #             CD <- f_DragCoefficient(DragMethod, Re, Psi, CSF, kS, kN)
  #             v_s_new <- sqrt(4 / 3 * d_eq / CD * ((rhoParticle - rhoFluid) / rhoFluid) * GN) * Cunningham
  #             RSS <- (v_s - v_s_new) ^ 2
  #             return(RSS)}
  #           result <- optimize(RSS_function, interval = c(0, 10), tol = 1e-10)
  # 
  #           return(result$minimum)
  #           },
  #           NA
  # )

}
