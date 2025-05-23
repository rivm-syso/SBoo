% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_Ksw.R
\name{f_Ksw}
\alias{f_Ksw}
\title{soil water partitioning coefficient}
\usage{
f_Ksw(Kow, pKa, CorgStandard, a, b, ChemClass, RHOsolid, alt_form, Ksw_orig)
}
\arguments{
\item{Kow}{octanol water partitioning coefficient [-]}

\item{CorgStandard}{Standard mass FRACTION organic carbon in soil/sediment [-]}

\item{a}{see QSAR table}

\item{b}{see QSAR table}

\item{ChemClass}{Class of chemical, see QSAR table (REACH, 2012)}

\item{RHOsolid}{density of "solids" [kg/m3]}

\item{alt_form}{the alternate form.}

\item{Ksw_orig}{soil water partitioning coefficient as present in the data [-]}
}
\value{
Ksw
}
\description{
Calculates the Ksw based on QSARS (Franco & Trapp, 2008) (https://doi.org/10.1897/07-583.1)
}
