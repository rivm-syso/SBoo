#' @title Correction factor depth-dependency soil concentrations
#' @name f_CORRsoil
#' @description 
#' Correction factor depth dependent soil concentrations. 
#' This is used for correcting the rate of runoff and volatilisation 
#' due to the non-homogeneaus distribution across the top soil layer.
#' The default relevant depth is 0 cm and the soil penetration depth is 10 cm.
#' The approach to derive this correction factor can be obtained from Hollander et al. (2004).
#' @param vertDistance depth of regional, continental and global soils [m]
#' @param penetration_depth_s (get discription from Hollander et al. (2004))
#' @param relevant_depth_s (get discription from Hollander et al. (2004))
#' @return Correction factor depth dependent soil concentration
#' @export
f_CORRsoil <- function (VertDistance, relevant_depth_s, penetration_depth_s){

    return(exp((-1/penetration_depth_s)*relevant_depth_s)*(1/penetration_depth_s) * VertDistance/(1-exp((-1/penetration_depth_s) * VertDistance))) #[-])

  }
