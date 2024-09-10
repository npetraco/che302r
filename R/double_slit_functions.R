#--------------------------------------------
#' @title Double slit intensity function
#' @description This is a function to evaluate the double slit intensity function
#'
#' @param I.max The largest wavenumber to use.
#'
#' @details This is a function to evaluate the double slit intensity function
#'
#' @return XXXX
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
