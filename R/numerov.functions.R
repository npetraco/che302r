#--------------------------------------------
# Numerov's recursive formula for the wave function
# (f is the wave function. np1/nm1 means n plus/minus one)
#--------------------------------------------
f.point.plus.one<-function(f.point,f.point.minus.one,increment,G.point.plus.one,G.point,G.point.minus.one) {
  return((2*f.point - f.point.minus.one + (5*G.point*f.point*increment^2)/6 + (G.point.minus.one*f.point.minus.one*increment^2)/12)/(1 - (G.point.plus.one*increment^2)/12))
}


#--------------------------------------------
#G function:
#--------------------------------------------
G<-function(disc.V, an.Energy) {
  return(2*(disc.V - an.Energy))
}


#--------------------------------------------
#' @title Potential energy function for use with Numerov's procedure
#' @description Potential energy function
#'
#' @param xax The x-axis over which to find a solution to the Schrodinger equation with the requested potential
#' @param potential.name Name of the potential energy function to insert into the Schrodinger equation. Choices are: "box", "harmonic", "anharmonic", "radial"
#' @param l The l quantum number if the "radial" potential energy function is requested
#'
#' @details Potential energy functions for use with Numerov's procedure to solve Schrodinger
#' equations implemented in this package. If Laguerre solutions are desired for orbital
#' calculations, choose the "radial" potential energy function and specify an l quantum number.
#'
#'
#' @return Values of the potential along the specified x-axis.
#'
#' @references Quantum Chemistry by Ira Levine. Any edition.
#'
#' @examples
#' # Potential energy function for particle in a box:
#' # Domain (x-axis):
#' dx    <- 0.01                            # Increment on x-axis
#' x.min <- (0)                             # Start of x-axis. Well into (left) classically forbidden region
#' x.max <- 1.5                             # End of x-axis. Well into the (right) classically forbidden region
#' x     <- seq(from=x.min,to=x.max,by=dx)  # x-axis
#'
#' potential<-"box"
#' plot(x,V(x,potential))
#'
#' # Potential energy function for anharmonic oscillator:
#' #Domain (x-axis):
#' dx    <- 0.01
#' x.min <- (-0.8)
#' x.max <- 5
#' x     <- seq(from=x.min,to=x.max,by=dx)
#' plot(x,V(x,"anharmonic"))
#'
#--------------------------------------------
V<-function(xax,potential.name, l=NULL) {

  if( !(potential.name %in% c("box", "harmonic", "anharmonic", "radial")) ){
    #print("Error!")
    #print("Your choices are: box, harmonic, anharmonic or laguerre.")
    stop("Bad choice in PE function. Your choices are: box, harmonic, anharmonic or radial.")
  }

  if(potential.name=="box"){
    return(rep(0,length(xax)))         #Particle in a box V
  }

  if(potential.name=="harmonic"){
    return((1/2) * xax^2)            #Harmonic oscillator V
  }

  if(potential.name=="anharmonic"){
    return(151.29 * (1-exp(-xax))^2) #Morse V, approx anharmonic oscillator
    #return(7.61 * (1-exp(-(xax-74.1) * 0.0193))^2)
  }

  if(potential.name=="radial"){
    return(   0.5*((l*(l+1))/xax^2) -(1/xax))  #... used to send in l
  }

}


#--------------------------------------------
#' @title Numerov's procedure
#' @description Numerov's procedure to solve simple 1D Schodinger equations
#'
#' @param xaxis The x-axis over which to find a solution to the Schrodinger equation with the requested potential
#' @param dxaxis The "differential" increment for the x-axis
#' @param PE.function.name Potential energy function name. Choices are: "box", "harmonic", "anharmonic", "radial"
#' @param nodes.in.state Number of nodes that should be in the wave function.
#' @param E.guess Guess for the energy of the wave function
#' @param num.iterations Number of iterations to run the procedure
#' @param delay.time Set this if you want to slow down the procedure and watch the solution evolve over the iterations.
#'
#' @details Numerov's procedure to solve simple 1D Schrodinger equations. Works for any of
#' the above stated potential energy functions. See Ira Levine's Quantum Chemistry for a very
#' nice explanation of the algorithm.
#'
#'
#' @return A two column matrix. The first column is the x-axis over which the wave function
#' was obtained. The second column is the un-normalized wave function solution.
#'
#' @references Quantum Chemistry by Ira Levine. Any edition.
#'
#' @examples
#' # Solve the Schrodinger equation for particle in a box:
#' #Domain (x-axis):
#' dx    <- 0.01                           # Increment on x-axis
#' x.min <- (0)                            # Start of x-axis. Well into (left) classically forbidden region
#' x.max <- 1.5                            # End of x-axis. Well into the (right) classically forbidden region
#' x     <- seq(from=x.min,to=x.max,by=dx) # x-axis
#'
#' #Potential energy function
#' potential<-"box"
#' plot(x,V(x,potential))
#'
#' #Numerov's procedure:
#' state        <- 0                # State you want, starting from 0. It is an integer: 0,1,2,3,...etc. I.E., number of nodes
#' Guess.Energy <- (0)              # Initial guess for energy of the state
#' max.iter     <- 40
#' psi.info     <- numerov.procedure(x,dx,potential,state,Guess.Energy,max.iter,delay.time=0.00)
#'
#' #Approximately Normalize Psi:
#' Npsi.info <- approx.normalize(psi.info)
#' plot(Npsi.info[,1], Npsi.info[,2], xlab="x", ylab="N*psi(x)", typ="l")
#'
#'
#' # Solve the Schrodinger equation for radial wave functions (Basically Laguerre functions):
#' #Domain (x-axis):
#' dr    <- 0.1                            # Increment on r-axis
#' r.min <- 1e-15                          # Start of r-axis.
#' r.max <- 180                            # End of r-axis.
#' r     <- seq(from=r.min,to=r.max,by=dr) # r-axis
#'
#' #Numerov's procedure to solve the (almost) Radial Schrodinger Equation:
#' n            <- 6      # n = 1,2,3,...
#' l            <- 4      # l = 0, ... , n-1 (e.g. s,p,d,f,g,.....)
#' Guess.Energy <- (-0.5) # Initial guess for energy of the state
#' max.iter     <- 50
#' state        <- (n-l-1)
#' Fr.info      <- numerov.procedure(r,dr,"radial",state,Guess.Energy,max.iter,0.00)
#'
#' #Transform F(r) back into R(r), the proper radial wave function
#' psi.info     <- Fr.info
#' psi.info[,2] <- (r^-1)*Fr.info[,2]
#' plot(r, r^2 * psi.info[,2],typ="l")
#'
#' #Approximately Normalize Psi:
#' Npsi.info <- approx.normalize(psi.info, include.jacobian = TRUE)
#'
#' #Plot the normalized wave function:
#' plot(r, r^2*Npsi.info[,2], xlim=c(min(r),max(r)), ylab="N psi(x)", typ="l")
#'
#--------------------------------------------
numerov.procedure <- function(xaxis, dxaxis, PE.function.name, nodes.in.state, E.guess,num.iterations, delay.time=0.00){

  Elow<-NULL  #Initialization
  Ehigh<-NULL
  Er<-E.guess #First iteration Enenrgy

  for(iter in 1:num.iterations) {

    if(PE.function.name=="radial"){
      Gr<-G(V(xaxis,PE.function.name,l), Er)
    } else{
      Gr<-G(V(xaxis,PE.function.name), Er)
    }


    psir<-numeric(length(xaxis))  #Initialize a psir for the iteration
    psir[1]<-0
    psir[2]<-(1*10^(-4))

    #Compute the psir for this iteration
    num.nodes<-0
    for(i in 3:length(xaxis)) {
      psir[i] <- f.point.plus.one(psir[i-1], psir[i-2], dxaxis, Gr[i], Gr[i-1], Gr[i-2])
      if(psir[i-1]*psir[i]<0) {
        num.nodes<-num.nodes+1
      }
      #print(psir[i])
    }

    #Check the energy and update it:
    if(num.nodes<=nodes.in.state) {
      Elow<-c(Elow,Er)
      if(is.null(Ehigh)) {
        print("E too low, but no Ehigh. Add 1 to E.")
        Er<-(Er+1)
      } else {
        print("E too low. Find avg. of max Elow and min Ehigh")
        print(paste("E bracket: [",max(Elow),",",min(Ehigh),"]"))
        Er<-(max(Elow)+min(Ehigh))/2
      }

    } else {
      print("E too high. Average with max Elow")
      Ehigh<-c(Ehigh,Er)
      print(paste("E bracket: [",max(Elow),",",min(Ehigh),"]"))
      Er<-(max(Elow) + min(Ehigh))/2
    }

    print(paste("Iteration: ",iter,"   E: ",Er,"   # of nodes: ",num.nodes,sep=""))
    plot(xaxis,psir)
    Sys.sleep(delay.time)
  }
  return(cbind(xaxis,psir))
}


#--------------------------------------------
#' @title Approximate 1D normalization
#' @description Approximately normalize a 1D-wave function
#'
#' @param wf.info Output from numerov.procedure
#' @param include.jacobian For use with radial wave functions. To normalize these wave functions the Jacobian must be included (TRUE)
#' @param plotQ Whether or not to plot the un-normalized probability density. For debugging.
#'
#' @details Approximately normalize a 1D-wave function spit out of numerov.procedure.
#'
#'
#' @return A two column matrix. The first column is the x-axis over which the wave function
#' was obtained. The second column is the normalized wave function solution.
#'
#' @references Quantum Chemistry by Ira Levine. Any edition.
#'
#' @examples
#' cf. examples in numerov.procedure
#'
#--------------------------------------------
approx.normalize<-function(wf.info, include.jacobian=FALSE, plotQ=FALSE) {
  xx<-wf.info[,1]
  wfx<-wf.info[,2]

  #Fit a spline function to the wave function^2 values:
  if(include.jacobian == TRUE) {
    den.func<-splinefun(xx, 4*pi * xx^2 * wfx^2)
  } else {
    den.func<-splinefun(xx,wfx^2)
  }

  #Compute the integral needed in the normalization constant = Int Psi* Psi dx:
  intgpsp<-integrate(den.func,lower=min(xx), upper=max(xx))$value

  #Compute the normalization constant:
  N.const<-1/sqrt(intgpsp)
  if(plotQ==TRUE) {
    plot(xx,(den.func(xx))^2,ylab="|psi(x)|^2")
    print(paste("Approx. normalization const.:",N.const))
  }

  #Scale the wave function values with the normalization constant
  wf.info.normalized<-wf.info
  wf.info.normalized[,2]<-(N.const*wf.info.normalized[,2])

  return(wf.info.normalized)
}
