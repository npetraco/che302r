# library(rgl)
# library(orthopolynom)
# library(parallel)
# library(doMC)
# library(doSNOW)

#--------------------------------------------
#Convert spherical coordinates to cartesian coordinates
#--------------------------------------------
spherical.to.cartesian<-function(rl,tha,pha) {
  x<- rl * sin(tha) * cos(pha)
  y<- rl * sin(tha) * sin(pha)
  z<- rl * cos(tha)

  return(c(x,y,z))
}

#--------------------------------------------
#Convert cartesian coordinates to spherical coordinates
#--------------------------------------------
cartesian.to.spherical<-function(xx.in,yy.in,zz.in) {

  rl  <- sqrt(xx.in^2 + yy.in^2 + zz.in^2)
  tha <- acos(zz.in/rl)
  pha <- atan(xx.in/yy.in)

  return(c(rl,tha,pha))
}

#--------------------------------------------
#Sample the wave function's density and plot corresponding r,theta, phi (x,y,z) values
#Recipe adapted from D. Cromer, J Chem Educ. 45(10) 626-631 1968
#--------------------------------------------
sample.density2<-function(nqn, lqn, mqn, spherical.grid, num.samples) {

  r<-spherical.grid[,1]
  theta<-spherical.grid[,2]
  phi<-spherical.grid[,3]

  #den.vals<-Re(psi.func(r,theta,phi)^2) #The Imaginary parts should be 0, but still need to extract Re
  #den.vals<-orbital_density(nqn, lqn, mqn, spherical.grid) #Use the boost functions
  den.vals <- r^2 * sin(theta) * Re((Radial(n, l, r) * LegendreP(l, m, theta) * PhiF(m, phi))^2)
  den.max<-max(den.vals)                  #Approx max of the density
  pratio.vals<-den.vals/den.max           #Ratio for monte carlo selection

  count<-0
  num.iter.keep<-num.samples #For now and speed sample the number of samples. This is a big change from Cromers paper. Much faster but may be wrong.
  #Prune the grid around the most probable values:
  keep.grid<-NULL
  while(count<num.samples) {
    q.iter<-runif(length(pratio.vals))
    idxs<-which(pratio.vals>=q.iter)
    if(length(idxs>0)) {
      idx.pick<-sample(idxs,num.iter.keep)
      keep.grid<-rbind(keep.grid,spherical.grid[idx.pick,])
      count<-(count+length(idx.pick))
      #print(count)
    }
  }
  orbital<-t(sapply(1:num.samples,function(x){spherical.to.cartesian(keep.grid[x,1],keep.grid[x,2],keep.grid[x,3])}))

  return(orbital)

}

#--------------------------------------------
#
#Radial wave functions using the orthopoly package
#
#--------------------------------------------
#Generates the symbolic expression for the desired Laguerre polynomial
#--------------------------------------------
laguerreLgen<-function(nn,ll){

  poly.nom <- glaguerre.polynomials(nn, ll, normalized=F)
  poly.nom <- list(poly.nom[[length(poly.nom)]]) #to gen polynomial values orthopolynom func expexts a list

  return(poly.nom)

}

#--------------------------------------------
#Generate a numerical representation of the desired radial wavefunction
#--------------------------------------------
Radial<-function(nn,ll,rr){

  lag.poly <- laguerreLgen(nn-ll-1, 2*ll+1)
  #print(lag.poly)

  #Radial.vals <- exp(-rr/nn) * rr^ll * polynomial.values(lag.poly,rr)[[1]]
  Radial.vals <- sqrt(factorial(nn-ll-1)/factorial(nn+ll)) * exp(-rr/nn) * (2*rr/nn)^ll * (2/(nn)^2) * polynomial.values(lag.poly,2*rr/nn)[[1]]

  return(Radial.vals)

}

#--------------------------------------------
#
#Angular wave functions using the Numerical Recipies routine
#
#There is no associated Legendre polynomial or spherical harmonic in orthopolynom, so we have to do it ourselves
#--------------------------------------------
#Generates the symbolic expression for a Legendre polynomial
#--------------------------------------------
legendrePgen<-function(ll){

  poly.nom <- legendre.polynomials(ll, normalized=F)
  poly.nom <- list(poly.nom[[length(poly.nom)]]) #to gen polynomial values orthopolynom func expects a list

  return(poly.nom)

}

#--------------------------------------------
#Generates the symbolic expression for |m|^th derivative Legendre polynomial. Needed for Associated Legendre polynomials
#--------------------------------------------
mderivlegendrePgen<-function(ll,mm){

  polyd <- legendrePgen(l)
  if(abs(mm)>0){
    for(i in 1:abs(mm)){
      polyd <- polynomial.derivatives(polyd)
    }
  }

  return(polyd)
}
#--------------------------------------------
#A pure function for the associated Legendre polynomial
#--------------------------------------------
LegendrePfunc<-function(ll,mm){

  deriv.func <- polynomial.functions(mderivlegendrePgen(ll,mm))[[1]]

  Plmx <- function(xx){
    (-1)^abs(mm) * (1-xx^2)^(abs(mm)/2) * deriv.func(xx)
  }

  return(Plmx)

}

#--------------------------------------------
#Compute a numerical representation for the associated Legendre polynomial
#--------------------------------------------
LegendreP<-function(ll,mm,thetaa){

  xxx<-cos(thetaa)
  vals <- LegendrePfunc(ll,mm)(xxx)

  if(mm<0){
    vals <- (-1)^mm * factorial(ll-mm)/factorial(ll+mm) * vals
  }

  normLP <- sqrt(((2*ll + 1)*factorial(ll - abs(mm)))/(2*factorial(ll + abs(mm))))

  return(normLP * vals)

}

#--------------------------------------------
#Compute the (numerical) Phi wave function for the equitorial angle
#--------------------------------------------
PhiF<-function(mm,phii){

  return(1/sqrt(2*pi) * exp((mm*complex(real = 0, imaginary = 1))*phii))

}



#---------------------------------------------------
# MCMC subroutines
#---------------------------------------------------
#
#--------------------------------------------
# Generate spline versions of the squared orbital component functions:
# This will save us from having to re-evaluate using the slower polynomial functions during MCMC
#--------------------------------------------
splined.sqorb.functs <- function(n, l, m, r.max = 20, num.knots = 1000){

  r     <- seq(from=1e-6, to=r.max, length.out = num.knots)
  theta <- seq(from=0,    to=pi,    length.out = num.knots)
  phi   <- seq(from=0,    to=2*pi,  length.out = num.knots)

  R.sq.vals <- sapply(r, function(r){Radial(n, l, r)^2})
  R.sq.func <- splinefun(r, R.sq.vals)
  #plot(r, r^2*R.den.func(r),typ="l")      # Check: is tail of Radial part long enough?
  #plot(r, log(r^2*R.den.func(r)),typ="l") # Check: is number of Radial nodes correct?

  Th.sq.vals <- sapply(theta, function(theta){sin(theta) * LegendreP(l, m, theta)^2})
  Th.sq.func <- splinefun(theta, Th.sq.vals)
  #plot(theta, sin(theta)*Th.den.func(theta), typ="l")
  #plot(theta, log(sin(theta)*Th.den.func(theta)), typ="l")

  Ph.sq.vals <- sapply(phi, function(phi){Re(PhiF(m, phi))^2})
  Ph.sq.func <- splinefun(phi, Ph.sq.vals)
  #plot(phi, Ph.den.func(phi), typ="l")
  #plot(phi, log(Ph.den.func(phi)), typ="l")

  spline.orb.func.info <- list(
    r,
    theta,
    phi,
    R.sq.func,
    Th.sq.func,
    Ph.sq.func
  )

  names(spline.orb.func.info) <- c(
    "r.axis",
    "theta.axis",
    "phi.axis",
    "R.sq.func",
    "Theta.sq.func",
    "Phi.sq.func"
  )

  return(spline.orb.func.info)

}

#--------------------------------------------
# Next step proposal for Metropolis algorithm, generalized to any number of parameters.
# CAUTION: uses Gaussian proposal
#--------------------------------------------
proposal <- function(theta.curr, proposal.wid){sapply(1:length(theta.curr), function(xx){rnorm(n = 1, mean = theta.curr[xx], sd = proposal.wid[xx])})}

#--------------------------------------------
# Orbital Ansatz using spline functions
# Log Density ansatz. Pass in splined function info to hide from user. ****Looks yucky though!
#--------------------------------------------
ansatz <- function(a.theta, R.den.func, Th.den.func, Ph.den.func, rax.ext){

  # abs the density functions incase they go slightly negative

  if((a.theta[1] > 0) & (a.theta[1] <= rax.ext)){
    val1 <- log( a.theta[1]^2 * abs(R.den.func(a.theta[1]))  )  # Replaces: Radial(n, l, a.theta[1])^2
  } else {
    val1 <- log(0)
  }

  if( (a.theta[2] >= 0) & (a.theta[2] <= pi) ) {
    val2 <- log( sin(a.theta[2]) * abs(Th.den.func(a.theta[2])) ) # Replaces: (LegendreP(l, m, a.theta[2]))^2
  } else {
    val2 <- log(0)
  }

  if( (a.theta[3] >= 0) & (a.theta[3] <= 2*pi)){
    val3 <- log( abs(Ph.den.func(a.theta[3])) )                   # Replaces: Re(PhiF(m, a.theta[3]))^2
  } else {
    val3 <- log(0)
  }

  val <- val1 + val2 + val3

  return(val)
}

#--------------------------------------------
# Log Metropolis ratio. Pass in splined function info to hide code...... Looks yucky
#--------------------------------------------
r.metrop <- function(theta.curr, theta.prop, rsf, thsf, phsf, rext){
  val <- ansatz(theta.prop, rsf, thsf, phsf, rext) - ansatz(theta.curr, rsf, thsf, phsf, rext)
  if(is.nan(val)){
    val <- (-Inf)
  }
  return(val)
}

#--------------------------------------------
# MCMC Metropolis routine:
# Should execute for any number of dimensions for a parameter vector.
# Requires that an ansatz (a not necessarily normalized pdf) be defined
#--------------------------------------------
sampler.gen <- function(num.iter=100, theta.init=c(0.5,0.5), proposal.width=c(0.5,0.5), spline.func.info){

  theta.current     <- theta.init
  ansatz.sample     <- array(NA,c(num.iter, length(theta.init))) # to theta sample from ansatz
  ansatz.sample[1,] <- theta.current
  accept.ps         <- array(NA,num.iter-1) # capture all acceptance probs just for checks
  accept.indcs      <- array(NA,num.iter-1) # capture all acceptance indicators just for checks

  # To pass into r.metrop and then to ansatz routine: YUK!!!!!!
  rsqf  <- spline.func.info$R.sq.func
  thsqf <- spline.func.info$Theta.sq.func
  phsqf <- spline.func.info$Phi.sq.func
  rmx   <- max(spline.func.info$r.axis)

  for(i in 2:num.iter){

    #print(i)

    # suggest new position
    theta.proposal <- proposal(theta.curr = theta.current, proposal.wid = proposal.width)

    # Accept proposal. Note log(ratio) and also pass in splined density parts. Can we do this better??
    accept.ratio <- r.metrop(
      theta.curr = theta.current,
      theta.prop = theta.proposal, rsqf, thsqf, phsqf, rmx)

    acceptQ           <- log(runif(1)) < accept.ratio
    accept.ps[i-1]    <- min(1,exp(accept.ratio))     # Just for checking
    accept.indcs[i-1] <- acceptQ                      # Just for checking

    if(acceptQ) {
      # Update position
      theta.current <- theta.proposal
    }

    ansatz.sample[i,] <- theta.current

  }

  return(list(ansatz.sample, accept.ps, accept.indcs))

}

#--------------------------------------------
#Sample the wave function's density and plot corresponding r,theta, phi (x,y,z) values
#Now uses MCMC Metropolis algorithm.
#
# Wrapper to run MCMC for the orbital densities:
#--------------------------------------------
sample.density3a<-function(nqn, lqn, mqn, orb.func.info, num.samples=5000, num.thin=1, num.burnin=0) {

  mpi  <- sampler.gen(
    theta.init       = runif(3, min = 0, max = 1),
    proposal.width   = c(1,1,1),
    num.iter         = num.samples,
    spline.func.info = orb.func.info)

  mp   <- mpi[[1]] # Chain
  ap   <- mpi[[2]] # Acceptance probs
  ai   <- mpi[[3]] # Acceptance indicators

  # Remove burn-in and Thin:
  keep.idxs <- seq(num.burnin+1,nrow(mp),num.thin)
  mp2  <- mp[keep.idxs, ]

  # Transform to x, y, z coords:
  samp.coords <- t(sapply(1:nrow(mp2),function(xx){spherical.to.cartesian(mp2[xx,1],mp2[xx,2],mp2[xx,3])}))

  colnames(samp.coords) <- c("x","y","z")

  print(paste("The final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}

#--------------------------------------------
# Wrapper to run MCMC for the orbital densities. This one we can input the proposal jump widths:
#--------------------------------------------
sample.density3b<-function(nqn, lqn, mqn, orb.func.info, num.samples=5000, num.thin=1, num.burnin=0, jump.widths=c(1,1,1), printQ=FALSE) {

  mpi  <- sampler.gen(
    theta.init       = runif(3, min = 0, max = 1),
    proposal.width   = jump.widths,
    num.iter         = num.samples,
    spline.func.info = orb.func.info)

  mp   <- mpi[[1]] # Chain
  ap   <- mpi[[2]] # Acceptance probs
  ai   <- mpi[[3]] # Acceptance indicators

  if(printQ==TRUE) {
    perc.1s    <- length(which(ap==1))/length(ap) * 100
    accpt.rate <- sum(ai)/length(ai) * 100 # Acceptance rate
    print(paste0("Higher move percentage: ", round(perc.1s,2), "%"))
    print(paste0("Acceptance rate:        ", round(accpt.rate,2), "%"))
  }

  # Remove burn-in and Thin:
  keep.idxs <- seq(num.burnin+1,nrow(mp),num.thin)
  mp2  <- mp[keep.idxs, ]

  # if(printQ==TRUE) {
  #   par(mfrow=c(3,1))
  #   acf(mp2[,1])
  #   acf(mp2[,2])
  #   acf(mp2[,3])
  # }

  # Transform to x, y, z coords:
  samp.coords <- t(sapply(1:nrow(mp2),function(xx){spherical.to.cartesian(mp2[xx,1],mp2[xx,2],mp2[xx,3])}))
  samp.coords <- cbind(samp.coords, mp2)

  colnames(samp.coords) <- c("x","y","z","r","theta","phi")

  print(paste("The final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}

#--------------------------------------------
# VERY stripped down wrapper (to a wrapper: sample.density3b) to run MCMC. CHECK THAT WE GOT ALL THE NODES THOUGH! MCMC by itself can easily miss
#--------------------------------------------
sample.orbital.density.mcmc<-function(nqn, lqn, mqn, orb.func.info, printQ=FALSE) {

  print("====== Normal pass ======")
  normal.pass.samp.size <- 187500

  samp.coords <- NULL

   for(i in 1:4) {

     print(paste("---- Starting Chain:", i, "----"))

     samp.coords <- rbind(samp.coords,
                          sample.density3b(nqn, lqn, mqn, orb.func.info,
                                           num.samples = normal.pass.samp.size,
                                           num.thin    = 150,
                                           jump.widths = c(2.38,2.38,2.38),
                                           num.burnin  = 100,
                                           printQ)
                          )
     print(paste("---- Chain:", i, "Done ----"))
  }

  colnames(samp.coords) <- c("x","y","z","r","theta","phi")

  print(paste("The total final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}

#--------------------------------------------
# VERY stripped down wrapper (to a wrapper: sample.density3b) to run MCMC. CHECK THAT WE GOT ALL THE NODES THOUGH! MCMC by itself can easily miss
# Simple parallel version: Run the chains in separate processes
#--------------------------------------------
sample.orbital.density.mcmc.parallel<-function(nqn, lqn, mqn, orb.func.info, printQ=FALSE, num.processes=1) {

  registerDoMC(num.processes)

  print(paste("Using",getDoParWorkers(), "processes."))

  print("====== Normal pass ======")
  normal.pass.samp.size <- 187500

  samp.coords <- foreach(i=1:4, .combine="rbind") %dopar% {
    sample.density3b(nqn, lqn, mqn, orb.func.info,
                     num.samples = normal.pass.samp.size,
                     num.thin    = 150,
                     jump.widths = c(2.38,2.38,2.38),
                     num.burnin  = 100,
                     printQ)
  }


  # x.loc <- density.samples2.loc[,1]
  # y.loc <- density.samples2.loc[,2]
  # z.loc <- density.samples2.loc[,3]
  #
  # samp.coords <- cbind(x.loc, y.loc, z.loc)
  #
  # colnames(samp.coords) <- c("x","y","z")

  colnames(samp.coords) <- c("x","y","z","r","theta","phi")

  print(paste("The total final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}


#--------------------------------------------
# Sample the wave function's density using AIS
# Recipe adapted from https://wiseodd.github.io/techblog/2017/12/23/annealed-importance-sampling/
#--------------------------------------------
sample.orbital.density.ais<-function(nqn, lqn, mqn, orb.func.info, printQ=FALSE) {

  rsqf  <- orb.func.info$R.sq.func
  thsqf <- orb.func.info$Theta.sq.func
  phsqf <- orb.func.info$Phi.sq.func
  rmx   <- max(orb.func.info$r.axis)

  # Work on the log-scale:

  # Target density f(0)
  logf0 <- function(an.x){

    loglik0 <-
      log(an.x[1]^2    * abs(rsqf(an.x[1])) )  +  # abs-ing incase splines went a little negative
      log(sin(an.x[2]) * abs(thsqf(an.x[2])) ) +
      log(               abs(phsqf(an.x[3])) )

    return(loglik0)
  }

  # Starting (Importance) density f(n)
  logfn <- function(an.x){

    return(
      dunif(an.x[1], min = 0, max = rmx,  log = T) +
        dunif(an.x[2], min = 0, max = pi,   log = T) +
        dunif(an.x[3], min = 0, max = 2*pi, log = T)
    )

  }

  # Generates an initial starting point for a sample from importance distribution, x ~ p(n)
  pn <- function(){

    return(
      c(runif(1, min = 0, max = rmx),
        runif(1, min = 0, max = pi),
        runif(1, min = 0, max = 2*pi))
    )

  }

  # Intermediate (annealing) densities
  logfj <- function(an.x, a.beta){
    return(a.beta * logf0(an.x) + (1-a.beta) * logfn(an.x))
  }

  # Transition kernel to bridge from starting density to target density with short MCMCa
  TransK <- function(an.x, a.beta, n.steps=10){
    for(t in 1:n.steps) {

      # Proposal
      x.prime <- c(NA,NA,NA)

      x.prime[1] <- an.x[1] + runif(1, min = -rmx, max = rmx)
      if(x.prime[1] < 0) {
        x.prime[1] <- x.prime[1] %% rmx
      } else if(x.prime[1] > rmx) {
        x.prime[1] <- x.prime[1] %% rmx
      }

      x.prime[2] <- an.x[2] + runif(1, min = -pi, max = pi)
      if(x.prime[2] < 0) {
        x.prime[2] <- x.prime[2] %% pi
      } else if(x.prime[2] > pi){
        x.prime[2] <- x.prime[2] %% pi
      }

      x.prime[3] <- an.x[3] + runif(1, min = -2*pi, max = 2*pi)
      if(x.prime[3] < 0) {
        x.prime[3] <- x.prime[3] %% 2*pi
      } else if(x.prime[3] > 2*pi){
        x.prime[3] <- x.prime[3] %% 2*pi
      }

      # log Acceptance prob
      log.accp <- logfj(x.prime, a.beta) - logfj(an.x, a.beta)

      if(is.na( log.accp )){
        print("Problem!!!!!!!")
        print(paste("beta:", a.beta))
        print(paste("x:", an.x))
        print(paste("log f_j(x, beta)", lf_j(an.x, a.beta)))
        print(paste("x.prime:", x.prime))
        print(paste("log f_j(x.prime, beta)", lf_j(x.prime, a.beta)))
      }

      out.x <- an.x
      if(log(runif(1)) < log.accp){
        out.x <- x.prime
      }

    }

    return(out.x)

  }

  # Do the sampling
  n.beta   <- 50 # Bridge over this many intermediate densities
  beta.seq <- seq(from = 0, to = 1, length.out = n.beta) # Inverse temp sequence

  n.samples   <- 2500 # Sample size sought from target density
  the.samples <- array(NA, c(n.samples,3))
  the.weights <- numeric(n.samples)

  # Do the sampling
  for(i in 1:n.samples){

    # Sample initial point from pn(x)
    xi <- pn()
    wi <- 1

    for(j in 2:length(beta.seq)){

      # Transition
      xi <- TransK(xi, beta.seq[j], n.steps=5)

      # Compute weight in log space (log-sum):
      wi <- wi + ( logfj(xi, beta.seq[j]) - logfj(xi, beta.seq[j-1]) )
    }

    if(printQ == TRUE) {
      print(paste("Done sample:", i))
    }

    the.samples[i,] <- xi
    the.weights[i] <- exp(wi)  # Transform back using exp

  }

  # Transform to x, y, z coords:
  samp.coords <- t(sapply(1:nrow(the.samples),function(xx){spherical.to.cartesian(the.samples[xx,1],the.samples[xx,2],the.samples[xx,3])}))
  samp.coords <- cbind(samp.coords, the.samples, the.weights)

  colnames(samp.coords) <- c("x","y","z","r","theta","phi", "wt")

  print(paste("The final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}

#--------------------------------------------
# Sample the wave function's density using AIS, but now do in parallel
#--------------------------------------------
sample.orbital.density.ais.parallel<-function(nqn, lqn, mqn, orb.func.info, printQ=FALSE, num.processes=1) {

  rsqf  <- orb.func.info$R.sq.func
  thsqf <- orb.func.info$Theta.sq.func
  phsqf <- orb.func.info$Phi.sq.func
  rmx   <- max(orb.func.info$r.axis)

  # Work on the log-scale:

  # Target density f(0)
  logf0 <- function(an.x){

    loglik0 <-
      log(an.x[1]^2    * abs(rsqf(an.x[1])) )  +  # abs-ing incase splines went a little negative
      log(sin(an.x[2]) * abs(thsqf(an.x[2])) ) +
      log(               abs(phsqf(an.x[3])) )

    return(loglik0)
  }

  # Starting (Importance) density f(n)
  logfn <- function(an.x){

    return(
      dunif(an.x[1], min = 0, max = rmx,  log = T) +
        dunif(an.x[2], min = 0, max = pi,   log = T) +
        dunif(an.x[3], min = 0, max = 2*pi, log = T)
    )

  }

  # Generates an initial starting point for a sample from importance distribution, x ~ p(n)
  pn <- function(){

    return(
      c(runif(1, min = 0, max = rmx),
        runif(1, min = 0, max = pi),
        runif(1, min = 0, max = 2*pi))
    )

  }

  # Intermediate (annealing) densities
  logfj <- function(an.x, a.beta){
    return(a.beta * logf0(an.x) + (1-a.beta) * logfn(an.x))
  }

  # Transition kernel to bridge from starting density to target density with short MCMCa
  TransK <- function(an.x, a.beta, n.steps=10){
    for(t in 1:n.steps) {

      # Proposal
      x.prime <- c(NA,NA,NA)

      x.prime[1] <- an.x[1] + runif(1, min = -rmx, max = rmx)
      if(x.prime[1] < 0) {
        x.prime[1] <- x.prime[1] %% rmx
      } else if(x.prime[1] > rmx) {
        x.prime[1] <- x.prime[1] %% rmx
      }

      x.prime[2] <- an.x[2] + runif(1, min = -pi, max = pi)
      if(x.prime[2] < 0) {
        x.prime[2] <- x.prime[2] %% pi
      } else if(x.prime[2] > pi){
        x.prime[2] <- x.prime[2] %% pi
      }

      x.prime[3] <- an.x[3] + runif(1, min = -2*pi, max = 2*pi)
      if(x.prime[3] < 0) {
        x.prime[3] <- x.prime[3] %% 2*pi
      } else if(x.prime[3] > 2*pi){
        x.prime[3] <- x.prime[3] %% 2*pi
      }

      # log Acceptance prob
      log.accp <- logfj(x.prime, a.beta) - logfj(an.x, a.beta)

      if(is.na( log.accp )){
        print("Problem!!!!!!!")
        print(paste("beta:", a.beta))
        print(paste("x:", an.x))
        print(paste("log f_j(x, beta)", lf_j(an.x, a.beta)))
        print(paste("x.prime:", x.prime))
        print(paste("log f_j(x.prime, beta)", lf_j(x.prime, a.beta)))
      }

      out.x <- an.x
      if(log(runif(1)) < log.accp){
        out.x <- x.prime
      }

    }

    return(out.x)

  }

  #-----------------
  # Do the sampling
  #-----------------
  n.beta    <- 50 # Bridge over this many intermediate densities
  beta.seq  <- seq(from = 0, to = 1, length.out = n.beta) # Inverse temp sequence
  n.samples <- 5000 # Sample size sought from target density

  # Sample function for parallel loop kernel. Every execution of this produces a sample:
  a.sample <- function() {

    # Sample initial point from pn(x)
    xi <- pn()
    wi <- 1

    for(jj in 2:length(beta.seq)){

      # Transition
      xi <- TransK(xi, beta.seq[jj], n.steps=5)

      # Compute weight in log space (log-sum):
      wi <- wi + ( logfj(xi, beta.seq[jj]) - logfj(xi, beta.seq[jj-1]) )
    }

    return(c(xi, exp(wi)))

  }

  # Run parallel loop here with progress bar
  cl <- makeSOCKcluster(num.processes)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=n.samples, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  the.samples <-
    foreach(i=1:n.samples,  .options.snow=opts, .combine='rbind') %dopar% {
      a.sample()
    }
  close(pb)
  stopCluster(cl)

  # Transform to x, y, z coords:
  samp.coords <- t(sapply(1:nrow(the.samples),function(xx){spherical.to.cartesian(the.samples[xx,1],the.samples[xx,2],the.samples[xx,3])}))
  samp.coords <- cbind(samp.coords, the.samples)

  colnames(samp.coords) <- c("x","y","z","r","theta","phi", "wt")

  print(paste("The final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}
