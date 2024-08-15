#' @title SB1Solve
#' @name SB1Solve
#' @description solve system of 1rst order k(i,j) and emissions,
#'  by solving v = 0
#' @param ParentModule SBcore
#' @param tol tolerance for accepting as steady state
#' @return States (i) (=mass)
#' export
SB1Solve = function(ParentModule, tol=1e-30) {
  
  SB.K = ParentModule$SB.k
  vEmis = ParentModule$emissions
  
  solve(SB.K, -vEmis, tol = tol)
}
