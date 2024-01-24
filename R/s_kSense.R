#' @title kSense
#' @name kSense
#' @description sensitivity Mass 4 kaas,
#' #logged sd / intervals ??
#'  by changing kaas what is the change in Mass
#' @param ParentModule SBCore
#' @param knamesdist description
#' @param samplesize description
#' @return States (i) (=mass)
#' @export
kSense = function(ParentModule, knames = NULL, tol=1e-30) {
  
  kaas <- ParentModule$myCore$kaas
  #check knames with kaas
  if (is.null(knames)) {
    knames <- unique(kaas$process)
  } else {
    stopifnot(all(knames %in% unique(kaas$process)))
  }
  # other preps; basic solution:
  SB.K = ParentModule$SB.k
  vEmis = ParentModule$emissions
  basicSolv <- solve(SB.K, -vEmis, tol = tol)
  
  aslist <- list()
  for (processNm in knames){
    kaas$k[kaas$process %in% knames] <- kaas$k[kaas$process %in% knames] * 1.01
    SB.K = ParentModule$PrepKaasM(kaas)
    aslist[[processNm]] <- solve(SB.K, -vEmis, tol = tol)
  }
  
  # rowsaslist <- list()
  # for (irow in 1:nrow(knamesdist)) {
  #   rowsaslist[rowsaslist] <- switch(knamesdist$distr,
  #     "rnorm", rnorm(samplesize, knamesdist$mean[irow], knamesdist$sd[irow]),
  #     "runif", runif(samplesize, knamesdist$Minm[irow], knamesdist$Maxm[irow]))
  # }
  # kaassamples <- docall(rbind, rowsaslist)
  
  do.call(rbind, aslist)
  
}
