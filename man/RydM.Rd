% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constants.R
\name{RydM}
\alias{RydM}
\title{"Corrected" Rydberg constant for any hydrogen like atom}
\usage{
RydM(at.mass, mass.unit = "amu", en.unit = "cm-1")
}
\arguments{
\item{at.mass}{atomic nucleus' mass}

\item{mass.unit}{"amu" or "kg"}

\item{en.unit}{Rydberg constant unit desired: "cm-1" (inverse centimeters), "m-1" (inverse meters) or J ("Joules)}
}
\value{
The corrected hydrogen-like atom Rydberg constant
}
\description{
"Corrected" Rydberg constant for any hydrogen like atom
}
\details{
The Rinf (or RinfJ) Rydberg constant assumes the center of mass for an nucleus-electron system is
directly at the nucleus (i.e. nucleus is clamped). To correct for this multiply by mu = me*mnuc/(me+mnuc)
and divide by me.
}
\examples{
RydM(mP, mass.unit = "kg", en.unit = "cm-1")
RH

RydM(1, mass.unit = "amu", en.unit = "cm-1")
RH # a tad different because 1 amu != mP or mN

RydM(2*mP+2*mN, mass.unit = "kg", en.unit = "cm-1") # Rydberg for Helium-4 +1 ion
}
