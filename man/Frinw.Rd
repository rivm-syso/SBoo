% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/v_FRinw.R
\name{Frinw}
\alias{Frinw}
\alias{FRinw}
\title{Fraction of chemical truly dissolved in water phase fresh water (relevant to Molecular species)}
\usage{
FRinw(
  FRorig_spw,
  FRACw,
  FRACa,
  FRACs,
  Kp,
  all.rhoMatrix,
  KpCOL,
  Kacompw,
  SUSP,
  COL,
  Matrix
)
}
\arguments{
\item{FRorig_spw}{Fraction original species in porewater of soil [-]}

\item{FRACw}{Fraction water in compartment [-]}

\item{FRACa}{Fraction air in compartment [-]}

\item{FRACs}{Fraction solids in compartment [-]}

\item{Kp}{soil/water partitioning coefficient [-]}

\item{Kacompw}{Dimensionless air/water partitioning coefficient of original species at compartment temperature [-]}

\item{SUSP}{Concentration suspended matter in water [mg.L-1]}

\item{COL}{Concentration of colloidal organic matter in water [mg.L-1]}

\item{Matrix}{Type of compartment matrix (soil, water, sediment or air)}

\item{rhoMatrix}{density of the matrix [kg.m-3]}
}
\description{
Everything that is not attached
}
