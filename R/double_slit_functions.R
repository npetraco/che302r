#--------------------------------------------
#' @title Double slit intensity function
#' @description This is a function to evaluate the double slit intensity function
#'
#' @param theta        direction angle between -pi/2 and pi/2 (-90 deg, 90deg)
#' @param d.slit.sep   slit separation
#' @param a.slit.width slit width
#' @param lambda       wavelength
#' @param I.max  largest intensity to use
#'
#' @details This is a function to evaluate the double slit intensity function.
#'
#' @return Double slit intensity.
#'
#'
#' @examples
#' library(che302r)
#'
#' dsi(pi/2, 10, 0.1, 650)
#'
#' #theta <- seq(from=-pi/2, to=pi/2, length.out=1000)
#' theta <- seq(from=-pi/4, to=pi/4, length.out=1000)
#' lam <- 650
#' d <- 6*lam
#' a <- 2*lam
#' I <- dsi(theta, d.slit.sep = d, a.slit.width = a, lambda = lam)
#' plot(theta*180/pi, I, typ="l", xlab="theta (deg)")
#'
#--------------------------------------------
dsi <- function(theta, d.slit.sep, a.slit.width, lambda, Imax=1) {

  a1    <- d.slit.sep*pi*sin(theta)/lambda
  a2    <- a.slit.width*pi*sin(theta)/lambda
  I.loc <- Imax * (cos(a1))^2 * (sin(a2)/a2)^2

  return(I.loc)
}


#--------------------------------------------
#' @title Single slit intensity function
#' @description This is a function to evaluate the single slit intensity function
#'
#' @param theta        direction angle between -pi/2 and pi/2 (-90 deg, 90deg)
#' @param a.slit.width slit width
#' @param lambda       wavelength
#' @param I.max  largest intensity to use
#'
#' @details This is a function to evaluate the single slit intensity function. It isn't really necessary
#' because we can just make the slit separation (d.slit.sep) 0 in dsi() to get ssi(). It's here for
#' convenience.
#'
#' @return Single slit intensity.
#'
#'
#' @examples
#' library(che302r)
#'
#' ssi(pi/2, 0.1, 650)
#'
#' theta <- seq(from=-pi/2, to=pi/2, length.out=1000)
#' #theta <- seq(from=-pi/4, to=pi/4, length.out=1000)
#' lam   <- 650
#' a     <- 3*lam
#' I     <- ssi(theta, a.slit.width = a, lambda = lam)
#' plot(theta*180/pi, I, typ="l", xlab="theta (deg)")
#'
#--------------------------------------------
ssi <- function(theta, a.slit.width, lambda, Imax=1) {

  a2    <- a.slit.width*pi*sin(theta)/lambda
  I.loc <- Imax * (sin(a2)/a2)^2

  return(I.loc)
}
