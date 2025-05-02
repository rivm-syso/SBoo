#' @title Settling Velocity Solver based on RSS
#' @name f_SetVelSolver
#' @description Calculates the settling velocity by minimizing  the residual sum of squares (RSS)
#' @param CD Drag Coefficient of a particle [-]
#' @param DragMethod Method of calculating the Drag Coefficient
#' @param Psi Shape factor, circularity/sphericity [-]
#' @param Re Reynolds number, as returned by the solver [-]
#' @param CSF Corey Shape Factor [-]
#' @param d_eq Equivalent spherical diameter of the particle [-]
#' @param DynViscFluidStandard Dynamic viscosity of liquid  (fraction of) compartment [kg.m-1.s-1]
#' @param rhoParticle Density of nanoparticle [kg.m-3]
#' @param rhoFluid Density of the water or air [kg m-3]
#' @param rad_species radius of the species [m]
#' @param GN gravitational force constant [m2 s-1]
#' @param Cunningham Cunningham coefficient, see f_Cunningham [-]
#' @return settling velocity [m/s] 
#' @export
#' 
f_SetVelSolver <- function(d_eq, Psi, DynViscFluidStandard, rhoParticle, rhoFluid, DragMethod, CSF, Matrix, rad_species, kS, kN) {
  # Define the RSS function to be minimized
  GN <- constants::syms$gn
  switch (Matrix,
          "water" = {
  
            RSS_function <- function(v_s) {
              Re <- d_eq * v_s *rhoFluid / DynViscFluidStandard
              CD <- f_DragCoefficient (DragMethod, Re, Psi, CSF, kS, kN)
              v_s_new <- sqrt(4 / 3 * d_eq / CD * ((rhoParticle - rhoFluid) / rhoFluid) * GN)
              RSS <- (v_s - v_s_new) ^ 2
              return(RSS)}
            result <- optimize(RSS_function, interval = c(-1, 1), tol = 1e-9)
            # # Extract the optimal settling velocity
            # v_s_final <- result$minimum
            # 
            # # Now compute Re using the final v_s
            # Re_final <- d_eq * v_s_final * rhoFluid / DynViscFluidStandard
            # 
            # # Print it
            # print(Re_final)
            
            return(result$minimum)
            },
          "air" = {
            RSS_function <- function(v_s) {
              Re <- d_eq * v_s *rhoFluid / DynViscFluidStandard
              Cunningham <- f_Cunningham(rad_species)
              CD <- f_DragCoefficient (DragMethod, Re, Psi, CSF, kS, kN)
              v_s_new <- sqrt(4 / 3 * d_eq / CD * ((rhoParticle - rhoFluid) / rhoFluid) * GN) * Cunningham
              RSS <- (v_s - v_s_new) ^ 2
              return(RSS)}
            result <- optimize(RSS_function, interval = c(0, 100), tol = 1e-9)
            
            return(result$minimum)
            },
            NA
  )       

}
