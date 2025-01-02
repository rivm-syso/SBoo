#' @title SB1Solve
#' @name SB1Solve
#' @description solve system of 1rst order k(i,j) and emissions,
#'  by solving v = 0
#' @param ParentModule SBcore
#' @param tol tolerance for accepting as steady state
#' @return States (i) (=mass)
#' export
SB1Solve = function(SB.K, emissions, tol=1e-30) {
  
  solve(SB.K, -emissions, tol = tol)
}
