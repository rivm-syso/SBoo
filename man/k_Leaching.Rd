% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k_Leaching.R
\name{k_Leaching}
\alias{k_Leaching}
\title{Leaching of particles due to infiltration into deeper soil layers}
\usage{
k_Leaching(
  FRACinf,
  RAINrate,
  VertDistance,
  SpeciesName,
  penetration_depth_s,
  Kscompw
)
}
\arguments{
\item{FRACinf}{Fraction infiltration #[-]}

\item{RAINrate}{Average precipitation #[m/s]}

\item{VertDistance}{depth compartment #[m]}

\item{SpeciesName}{considered species}

\item{penetration_depth_s}{penetration depth for correction factor, based on Hollander et al. (2007) https://doi.org/10.1016/j.envpol.2006.09.018 [m]}

\item{Kscompw}{dimensionless partitioning coefficient soil-water [-]}
}
\value{
Leaching of aggregated (A) or free (S) enp species from natural soil #[s-1]
}
\description{
Leaching of particles due to infiltration into deeper soil layers
}
