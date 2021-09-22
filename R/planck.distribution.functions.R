h  <- 6.626068e-34  # Planck's const
cl <- 299792458     # Speed of light
kB <- 1.3806503e-23 # Boltzmann's const

#-------------------------------------------------------------------------------------
#This is a function to plot the Plank distribution as a function of wavenumber.
#
#ARGUEMENTS:
#nut.min: The smallest wavenumber to use. Use something at least little larger than 0.
#nut.max: The largest wavenumber to use.
#Temp: Temperature in Kelvin
#--------------------------------------------------------------------------------------
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
#Now USE the function planck.distribution.nutilde below:
#planck.distribution.nutilde(10000,3000000,1500)


#-------------------------------------------------------------------------------------
#This is a function to plot the Plank distribution as a function of wavelength.
#
#ARGUEMENTS:
#lam.min: The smallest wavelength to use. Use something at least little larger than 0.
#lam.max: The largest wavelength to use.
#Temp: Temperature in Kelvin
#--------------------------------------------------------------------------------------
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
#Now USE the function planck.distribution.nutilde below:
#planck.distribution.lambda(10*10^-9,10000*10^-9,1500)
#planck.distribution.lambda(10*10^-9,12000*10^-9,1500)


#-------------------------------------------------------------------------------------
#This is a function to plot the Plank distribution as a function of frequency.
#
#ARGUEMENTS:
#nu.min: The smallest frequency to use. Use something at least little larger than 0.
#nu.max: The largest frequency to use.
#Temp: Temperature in Kelvin
#--------------------------------------------------------------------------------------
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


#-------------------------------------------------------------------------------------
#This is a function to plot the Plank distribution as a function of wavenumber, wavelength or frequency.
#
#ARGUEMENTS:
#typ: The x-axis scale: wavenumber (m^-1), wavelength (m) or frequency (s)
#x.min: The smallest wavelength to use. Use something at least little larger than 0.
#x.max: The largest wavelength to use.
#Temp: Temperature in Kelvin
#--------------------------------------------------------------------------------------
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
