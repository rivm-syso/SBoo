% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/v_Kaers.R
\name{Kaers}
\alias{Kaers}
\title{Aerosol solids - air Partition coefficient}
\usage{
Kaers(
  Kaw25,
  Kow,
  Corg,
  RhoCOL,
  Matrix,
  Pvap25,
  MaxPvap,
  Sol25,
  MW,
  T25,
  ChemClass
)
}
\arguments{
\item{Kaw25}{Dimensionless gas/water partitioning coefficient of the original species at 25C [-]}

\item{Kow}{octanol water partitioning coefficient [-]}

\item{Corg}{mass fraction organic carbon in soil/sediment [-]}

\item{RhoCOL}{density of accumulation mode aerosol particle [kg m-3]}

\item{Matrix}{type of compartment considered}

\item{Pvap25}{vapor pressure of original species at 25 C [Pa]}

\item{MaxPvap}{maximum vapor pressure of original species at 25 C [Pa]}

\item{Sol25}{water solubility of original species at 25 C [mol.m-3]}

\item{MW}{molar weight of original species [kg mol-1]}

\item{ChemClass}{chemical class, see QSAR table}
}
\value{
Kaers
}
\description{
partitioning coefficient between aerosol solids and air [-]
}
