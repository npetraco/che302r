% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/planck.distribution.functions.R
\name{planck.distribution.nutilde}
\alias{planck.distribution.nutilde}
\title{Plank distribution as a function of wavenumber.}
\usage{
planck.distribution.nutilde(nut.min, nut.max, Temp, plotQ = FALSE)
}
\arguments{
\item{nut.min}{The smallest wavenumber to use. Use something at least little larger than 0.}

\item{nut.max}{The largest wavenumber to use.}

\item{Temp}{The temperature in Kelvin}

\item{plotQ}{Whether or not to plot the distribution}
}
\value{
A two column matrix. The first column is the nu-tilde-axis over which the Planck
distribution was computed. The second column is the Planck distribution values
}
\description{
This is a function to plot the Plank distribution as a function of wavenumber (m^-1).
}
\details{
This is a function to plot the Plank distribution as a function of wavenumber.
NOTE: The units of wavenumber for this function are m^-1 and NOT cm^-1.
}
\examples{
#USE the function planck.distribution.nutilde below:
planck.distribution.nutilde(10000,3000000,1500)

}
