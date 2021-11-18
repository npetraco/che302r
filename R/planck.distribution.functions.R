h  <- 6.626068e-34  # Planck's const
cl <- 299792458     # Speed of light
kB <- 1.3806503e-23 # Boltzmann's const
hb <- h/(2*pi)      # Reduced Planck's const

#--------------------------------------------
#' @title Plank distribution as a function of wavenumber.
#' @description This is a function to plot the Plank distribution as a function of wavenumber (m^-1).
#'
#' @param nut.min The smallest wavenumber to use. Use something at least little larger than 0.
#' @param nut.max The largest wavenumber to use.
#' @param Temp The temperature in Kelvin
#' @param plotQ Whether or not to plot the distribution
#'
#' @details This is a function to plot the Plank distribution as a function of wavenumber.
#' NOTE: The units of wavenumber for this function are m^-1 and NOT cm^-1.
#'
#' @return A two column matrix. The first column is the nu-tilde-axis over which the Planck
#' distribution was computed. The second column is the Planck distribution values
#'
#'
#' @examples
#' #USE the function planck.distribution.nutilde below:
#' planck.distribution.nutilde(10000,3000000,1500)
#'
#--------------------------------------------
planck.distribution.nutilde<-function(nut.min, nut.max, Temp, plotQ=FALSE) {

  #Make a nu-tilde (nut) axis. This is the x-axis.
  nut <- seq(from=nut.min, to=nut.max, length.out=2500)

  #Planck's dist as a function of nu-tilde (nut). This is the y-axis.
  rho <- 2*h*(cl^2)*(nut^3) * (1/(exp((h*cl*nut)/(kB*Temp))-1))

  if(plotQ==TRUE){
    #Make the plot.
    plot(nut, rho, typ="l", xlab="nu-tilde (m^-1)", ylab="Intensity", main="Planck's distribution")
  }

  plank.dist.info <- cbind(
    nut,
    rho
  )
  colnames(plank.dist.info) <- c("nu.tilde", "rho")

  return(plank.dist.info)

}


#--------------------------------------------
#' @title Plank distribution as a function of wavelength.
#' @description This is a function to plot the Plank distribution as a function of wavelength.
#'
#' @param lam.min The smallest wavelength to use. Use something at least little larger than 0.
#' @param lam.max The largest wavelength to use.
#' @param Temp The temperature in Kelvin
#' @param plotQ Whether or not to plot the distribution
#'
#' @details This is a function to plot the Plank distribution as a function of wavelength.
#'
#' @return A two column matrix. The first column is the lambda-axis over which the Planck
#' distribution was computed. The second column is the Planck distribution values
#'
#'
#' @examples
#' #Now USE the function planck.distribution.lambda below:
#' planck.distribution.lambda(10e-9,10000e-9,1500)
#' planck.distribution.lambda(10e-9,12000e-9,1500)
#'
#--------------------------------------------
planck.distribution.lambda<-function(lam.min, lam.max, Temp, plotQ=FALSE) {

  #Make a nu-tilde (nut) axis. This is the x-axis.
  lam <- seq(from=lam.min, to=lam.max, length.out=2500)

  #Planck's dist as a function of nu-tilde (nut). This is the y-axis.
  #rho <- 2*h*(cl^2)*(nut^3) * (1/(exp((h*cl*nut)/(kB*Temp))-1))
  rho <- (2*h*(cl^2)/(lam^5)) * (1/(exp((h*cl)/(lam*kB*Temp))-1))

  if(plotQ==TRUE){
    #Make the plot.
    plot(lam, rho, typ="l", xlab="lambda (m)", ylab="Intensity", main="Planck's distribution")
  }

  plank.dist.info <- cbind(
    lam,
    rho
  )
  colnames(plank.dist.info) <- c("lambda", "rho")

  return(plank.dist.info)

}


#--------------------------------------------
#' @title Plank distribution as a function of frequency.
#' @description This is a function to plot the Plank distribution as a function of frequency.
#'
#' @param nut.min The smallest frequency to use. Use something at least little larger than 0.
#' @param nut.max The largest frequency to use.
#' @param Temp The temperature in Kelvin
#' @param plotQ Whether or not to plot the distribution
#'
#' @details This is a function to plot the Plank distribution as a function of frequency.
#'
#' @return A two column matrix. The first column is the nu-axis over which the Planck
#' distribution was computed. The second column is the Planck distribution values
#'
#'
#' @examples
#' #USE the function planck.distribution.nu below:
#' planck.distribution.nu(10e10, 6e14, 1500, plotQ = T)
#'
#--------------------------------------------
planck.distribution.nu<-function(nu.min, nu.max, Temp, plotQ=FALSE) {

  #Make a nu-tilde (nut) axis. This is the x-axis.
  nu <- seq(from=nu.min, to=nu.max, length.out=2500)

  #Planck's dist as a function of nu. This is the y-axis.
  #rho <- 2*h*(cl^2)*(nut^3) * (1/(exp((h*cl*nut)/(kB*Temp))-1))
  #rho <- (2*h*(cl^2)/(lam^5)) * (1/(exp((h*cl)/(lam*kB*Temp))-1))
  rho <- (2*h*(nu^3)/(cl^2)) * (1/(exp((h*nu)/(kB*Temp))-1))

  if(plotQ==TRUE){
    #Make the plot.
    plot(nu, rho, typ="l", xlab="nu (s^-1)", ylab="Intensity", main="Planck's distribution")
  }

  plank.dist.info <- cbind(
    nu,
    rho
  )
  colnames(plank.dist.info) <- c("nu", "rho")

  return(plank.dist.info)

}


#--------------------------------------------
#' @title Plank distribution as a function of wavenumber, wavelength or frequency.
#' @description This is a function to plot the Plank distribution as a function of wavenumber, wavelength or frequency.
#' @param typ   Scale units for witch to compute the Planck distribution. Choices are: "wavenumber", "wavelength", "frequency"
#' @param x.min The smallest domain value to use. Use something at least little larger than 0
#' @param x.max The largest domain to use
#' @param Temp The temperature in Kelvin
#' @param plotQ Whether or not to plot the distribution
#'
#' @details This wrapper is a function to plot the Plank distribution as a function of wavenumber, wavelength or frequency.
#'
#' @return A two column matrix. The first column is the x-axis over which the Planck
#' distribution was computed (m^-1, m or s^-1). The second column is the Planck distribution values
#'
#'
#' @examples
#' #USE the function planck.distribution below:
#' planck.distribution(typ = "wavenumber", 10000, 3000000,  1500, plotQ = T)
#' planck.distribution(typ = "wavelength", 10e-9, 10000e-9, 1500, plotQ = T)
#' planck.distribution(typ = "frequency",  10e10, 6e14,     1500, plotQ = T)
#'
#--------------------------------------------
planck.distribution<-function(typ, x.min, x.max, Temp, plotQ=FALSE) {

  if(typ=="wavenumber"){

    planck.info <- planck.distribution.nutilde(nut.min = x.min, nut.max = x.max, Temp = Temp, plotQ = plotQ)
    return(planck.info)

  } else if(typ=="wavelength") {

    planck.info <- planck.distribution.lambda(lam.min = x.min, lam.max = x.max, Temp = Temp, plotQ = plotQ)
    return(planck.info)

  } else if(typ=="frequency") {

    planck.info <- planck.distribution.nu(nu.min = x.min, nu.max = x.max, Temp = Temp, plotQ = plotQ)
    return(planck.info)

  } else {
    stop("Specify the typ: wavenumber, wavelength or frequency")
  }

}
