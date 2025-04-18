% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/v_FRACw.R
\name{FRACw}
\alias{FRACw}
\title{FRACw}
\usage{
FRACw(subFRACa, subFRACw, subFRACs, Matrix)
}
\arguments{
\item{subFRACa}{subfraction of air in a non-air compartment [-]}

\item{subFRACw}{subfraction of water in a non-water compartment [-]}

\item{subFRACs}{subfraction of solids in a non-soil, non-sediment compartment [-]}

\item{Matrix}{type of compartment}
}
\value{
FRACw
}
\description{
Fraction of water in any matrix, 
either subFRACw or, when the main matrix, remainder after substracting subFRACs + subFRACa
}
