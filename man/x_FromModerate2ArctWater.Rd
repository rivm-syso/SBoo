% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/x_FromModerate2ArctWater.R
\name{x_FromModerate2ArctWater}
\alias{x_FromModerate2ArctWater}
\title{Mixing of upper and deeper sea water layers (w2 and w3)}
\usage{
x_FromModerate2ArctWater(OceanCurrent, SubCompartName, ScaleName)
}
\arguments{
\item{OceanCurrent}{current implemented in the ocean [m3.s-1]}

\item{SubCompartName}{name of the subcompartment of the box at hand}

\item{ScaleName}{name of the scale of the box at hand}
}
\value{
Flow [m3.s-1]
}
\description{
Moderate scale waters flow into Arctic scale through surface ocean currents based on a set ocean current.
}
