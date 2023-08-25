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
kSense = function(ParentModule, knamesdist, samplesize) {
  
  ImplemDistr <- c("rnorm", "runif")
  #check (knamesdist for) anomilies
  stopifnot(all(c("kname", "distr", "sd") %in% names(knamesdist)))
  stopifnot(knamesdist$distr %in% ImplemDistr)
  kaas <- ParentModule$myCore$kaas
  #add k.Abbr if not present
  if (!"k.Abbr" %in% names(kaas)) {
    kaas$k.Abbr <- paste("k.", kaas$fromAbbr)
    NonDiag <- kaas$toAbbr != kaas$fromAbbr
    kaas$k.Abbr[NonDiag] <- paste(kaas$k.Abbr[NonDiag], kaas$toAbbr[NonDiag])
  }
  stopifnot(all(knamesdist$kname %in% kaas$Abbr))
  
  # other preps
  knamesInKaas <- match(knamesdist$k.Abbr, kaas$k.Abbr)
  kaasInKnames <- match(kaas$k.Abbr, knamesdist$k.Abbr)
  vEmis = ParentModule$PrepemisV()
  #Current value in kaas becomes mean in knamesdist
  # first and second parameter, depending on the distribution - type
  # mean, sd for rnorm; min, max for runif
  knamesdist$k.mean <- kaas$k[knamesInKaas]
  maxminmin = knamesdist$sd * sqrt(12)
  knamesdist$Minm <- knamesdist$k.mean - maxminmin / 2
  knamesdist$Maxm <- knamesdist$k.mean + maxminmin / 2
  rowsaslist <- list()
  for (irow in 1:nrow(knamesdist)) {
    rowsaslist[rowsaslist] <- switch(knamesdist$distr,
      "rnorm", rnorm(samplesize, knamesdist$mean[irow], knamesdist$sd[irow]),
      "runif", runif(samplesize, knamesdist$Minm[irow], knamesdist$Maxm[irow]))
  }
  kaassamples <- docall(rbind, rowsaslist)
  
  sresults <- list() #of frames to be rbind
  for (samplecol in 1:samplesize){
    kaas[kaasInKnames] <- kaassamples[,samplecol]
    SB.K = ParentModule$PrepKaasM(kaas)
    sresults[samplecol] <- solve(SB.K, -vEmis)
  }
  do.call(rbind, sresults)
  
}
