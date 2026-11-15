sample.orbital.density.mcmc_mod<-function(nqn, lqn, mqn, orb.func.info, pss=187500, thin=150, jumps=c(2.38,2.38,2.38), burnin=100, printQ=FALSE) {

  print("====== Normal pass ======")
  print(paste0("Thinned chain sample size will be: ", pss/thin))
  normal.pass.samp.size <- pss

  samp.coords <- NULL

  for(i in 1:4) {

    print(paste("---- Starting Chain:", i, "----"))

    samp.coords <- rbind(samp.coords,
                         sample.density3b(nqn, lqn, mqn, orb.func.info,
                                          num.samples = normal.pass.samp.size,
                                          num.thin    = thin,
                                          jump.widths = jumps,
                                          num.burnin  = burnin,
                                          printQ)
    )
    print(paste("---- Chain:", i, "Done ----"))
  }

  colnames(samp.coords) <- c("x","y","z","r","theta","phi")

  print(paste("The total final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}



sample.orbital.density.mcmc.parallel_mod<-function(nqn, lqn, mqn, orb.func.info, pss=187500, thin=150, jumps=c(2.38,2.38,2.38), burnin=100, printQ=FALSE, num.processes=1) {

  #registerDoMC(num.processes)
  #print(paste("Using",getDoParWorkers(), "processes."))

  print("====== Normal pass ======")
  print(paste0("Thinned chain sample sizes will be: ", (pss-burnin)/thin))
  print(paste0(num.processes, " chains will run for a total sample size = ", (pss - burnin)*num.processes/thin))
  normal.pass.samp.size <- pss

  # Run parallel loop here with progress bar
  cl <- makeSOCKcluster(num.processes)
  registerDoSNOW(cl)
  pb <- txtProgressBar(min=1, max=4, style=3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)

  samp.coords <-
    foreach(i=1:4, .export = c("sample.density3b"), .options.snow=opts, .combine='rbind') %dopar% {
      sample.density3b(nqn, lqn, mqn, orb.func.info,
                       num.samples = normal.pass.samp.size,
                       num.thin    = thin,
                       jump.widths = jumps,
                       num.burnin  = burnin,
                       printQ)
    }
  close(pb)
  stopCluster(cl)

  colnames(samp.coords) <- c("x","y","z","r","theta","phi")

  print(paste("The total final sample size is:", nrow(samp.coords)))

  return(samp.coords)

}


sample.orbital.density.ais_mod<-function(nqn, lqn, mqn, orb.func.info, num.samples=2500, num.beta=50, tkn.steps=5, printQ=FALSE) {

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
  TransK <- function(an.x, a.beta, n.steps){
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
  n.beta   <- num.beta # Bridge over this many intermediate densities
  beta.seq <- seq(from = 0, to = 1, length.out = n.beta) # Inverse temp sequence

  n.samples   <- num.samples # Sample size sought from target density
  the.samples <- array(NA, c(n.samples,3))
  the.weights <- numeric(n.samples)

  # Do the sampling
  for(i in 1:n.samples){

    # Sample initial point from pn(x)
    xi <- pn()
    wi <- 1

    for(j in 2:length(beta.seq)){

      # Transition
      xi <- TransK(xi, beta.seq[j], n.steps=tkn.steps)

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


sample.orbital.density.ais.parallel_mod<-function(nqn, lqn, mqn, orb.func.info, num.samples=5000, num.beta=50, tkn.steps=5, printQ=FALSE, num.processes=1) {

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
  TransK <- function(an.x, a.beta, n.steps){
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
  n.beta    <- num.beta # Bridge over this many intermediate densities
  beta.seq  <- seq(from = 0, to = 1, length.out = n.beta) # Inverse temp sequence
  n.samples <- num.samples # Sample size sought from target density

  # Sample function for parallel loop kernel. Every execution of this produces a sample:
  a.sample <- function() {

    # Sample initial point from pn(x)
    xi <- pn()
    wi <- 1

    for(jj in 2:length(beta.seq)){

      # Transition
      xi <- TransK(xi, beta.seq[jj], n.steps=tkn.steps)

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


sample.orbital.density_mod<-function(nqn, lqn, mqn, orb.func.info, printQ=FALSE, algorithm="ais", type="serial", mc.pss=187500, mc.thin=150, mc.jumps=c(2.38,2.38,2.38), mc.burnin=100, ais.num.samples=2500, ais.num.beta=50, ais.tkn.steps=5, num.chains=4, num.processes=1) {

  if(algorithm=="mcmc") {
    if(type=="serial") {
      samp.coords <- sample.orbital.density.mcmc(nqn = nqn, lqn = lqn, mqn = mqn, orb.func.info = orb.func.info, pss=mc.pss, thin=mc.thin, jumps=mc.jumps, burnin=mc.burnin, num.chains = num.chains, printQ = printQ)
    } else if(type=="parallel") {
      samp.coords <- sample.orbital.density.mcmc.parallel(nqn = nqn, lqn = lqn, mqn = mqn, orb.func.info = orb.func.info, pss=mc.pss, thin=mc.thin, jumps=mc.jumps, burnin=mc.burnin, num.chains=num.chains, printQ=FALSE, num.processes=num.processes)
    } else {
      stop("Choose serial or parallel")
    }
  } else if(algorithm=="ais") {
    if(type=="serial") {
      samp.coords <- sample.orbital.density.ais(nqn = nqn, lqn = lqn, mqn = mqn, orb.func.info = orb.func.info, printQ = printQ, num.samples=ais.num.samples, num.beta=ais.num.beta, tkn.steps=ais.tkn.steps)
    } else if(type=="parallel") {
      samp.coords <- sample.orbital.density.ais.parallel(nqn = nqn, lqn = lqn, mqn = mqn, orb.func.info = orb.func.info, printQ = printQ, num.samples=ais.num.samples, num.beta=ais.num.beta, tkn.steps=ais.tkn.steps, num.processes = num.processes)
    } else {
      stop("Choose serial or parallel")
    }
  } else {
    stop("Choose ais or mcmc")
  }

  return(samp.coords)

}
