#' @title Settling Velocity Solver based on RSS
#' @name f_SetVelSolver
#' @description Calculates the settling velocity by miDynViscWaterStandardmizing the residual sum of squares (RSS)
#' @param CD Drag Coefficient of a particle [-]
#' @param DragMethod Method of calculating the Drag Coefficient
#' @param Psi Shape factor, circularity/sphericity [-]
#' @param Re Reynolds number, as returned by the solver [-]
#' @param CSF Corey Shape Factor [-]
#' @return f_DragCoefficient 
#' @export
#' 
f_SetVelSolver <- function(d_eq, Psi, DynViscWaterStandard, rhoParticle, rhoWater, DragMethod, CSF) {
  # Define the RSS function to be miDynViscWaterStandardmized
  GN <- constants::syms$gn
  RSS_function <- function(v_s) {
    Re <- d_eq * v_s *rhoWater / DynViscWaterStandard
    CD <- f_DragCoefficient (DragMethod, Re, Psi, CSF)
    v_s_new <- sqrt(4 / 3 * d_eq / CD * ((rhoParticle - rhoWater) / rhoWater) * GN)
    RSS <- (v_s - v_s_new) ^ 2
    return(RSS)
  }
  
  # Use numerical optimizer to minimize the RSS function
  result <- optimize(RSS_function, interval = c(0, 1), tol = 1e-9)
  
  return(result$minimum)
}