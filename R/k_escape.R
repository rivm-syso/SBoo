#' @title escape
#' @name K_Escape
#' @description calculate k for escape from air compartment to stratosphere
#' @return k
#' @export
k_Escape <- function(t_half_Escape) 
  return(log(2)/(t_half_Escape*365*24*3600)) # no comment