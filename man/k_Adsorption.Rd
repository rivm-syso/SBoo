% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/k_Adsorption.R
\name{k_Adsorption}
\alias{k_Adsorption}
\title{Adsorption of chemical to across an interface}
\usage{
k_Adsorption(
  FRingas,
  FRinw,
  from.MTC_2sd,
  to.FRorig_spw,
  to.MTC_2w,
  from.MTC_2w,
  to.MTC_2a,
  from.MTC_2s,
  to.FRorig,
  Kacompw,
  to.Kscompw,
  to.Matrix,
  VertDistance,
  AreaLand,
  AreaSea,
  to.Area,
  all.FRorig,
  all.FRorig_spw,
  from.SubCompartName,
  to.SubCompartName,
  ScaleName,
  Test
)
}
\arguments{
\item{FRinw}{fraction of species in water,see FRinw()}

\item{Kacompw}{partitioning coefficient for air, see SBvariable}

\item{VertDistance}{Vertical distance of the compartment}

\item{AreaLand}{calculated land area of the scale considered, see SBvariable}

\item{AreaSea}{calculated sea area of the scale considered, see SBvariable}

\item{to.Area}{Area of the receiving compartent}

\item{ScaleName}{name of scale}

\item{Test}{Test = TRUE or FALSE, depending on if SBExcel version should be copied}

\item{MTC_2sd}{mass transfer coefficient to sediment, see MTC_2sd()}

\item{MTC_2w}{mass transfer coefficient to water, see SBvariable}

\item{MTC_2a}{mass transfer coefficient to air, see SBvariable}

\item{MTC_2s}{mass transfer coefficient to soil, see SBvariable}

\item{FRorig}{fraction of species in original form, see SBvariable}

\item{Kscompw}{partitioning coefficient for soil, see SB variable}

\item{Matrix}{Receiving compartment type}

\item{FRorig_spw}{fraction of original species in soil pore water}

\item{SubCompartName}{name of the subcompartname that is from/to}
}
\value{
The adsorption rate constant relevant for the receiving compartments soil, water or sediment [s-1]
}
\description{
The adsorption rate constant (k_Adsorption) is calculated based on the principles of gas absorption as described by Hollander
et al. (2004) 10.1080/10629360412331297470 & Hollander et al. (2007) https://doi.org/10.1016/j.envpol.2006.09.018
}
