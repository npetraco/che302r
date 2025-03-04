#--------------------------------------------
#' @title multinomial coefficient
#' @description This function computes the multinomial coefficient
#'
#' @param x         a macrostate, i.e. a vector of one-particle state occupancies
#' @param logQ      whether or not to return the log of the multinomial coefficient. Default is FALSE
#'
#' @details This is a function computes the multinomial coefficient.
#'
#' @return multinomial coefficient.
#'
#'
#' @examples
#' library(che302r)
#'
#' x <- c(3,2,8,0,0,1,0,0)
#' N <- sum(x)
#' factorial(N)/prod(factorial(x))
#' multinom(x)
#'
#' multinom(c(3,3,7,0,0,1,0,0))
#'
#' multinom(c(2,4,1,2,2,0,1,2), logQ = T)
#'
#--------------------------------------------
multinom <- function(x, logQ=F) {

  N   <- sum(x)
  val <- lfactorial(N) - sum(lfactorial(x))
  if(logQ == F) {
    val <- exp(val)
  }

  return(val)

}


#--------------------------------------------
#' @title Particle in a box energy
#' @description This function computes the particle in a box energy
#'
#' @param n vector of PIAB n quantum numbers: 1, 2, 3, .....
#' @param m particle mass in kg
#' @param L vector of box side lengths in meters
#'
#' @details This is a function computes the particle in a box energy. If n is a single positive integer
#' and L is a positive real number, the function computed the 1D particle in a box energy.
#'
#' @return a PIAB energy in joules.
#'
#'
#' @examples
#' library(che302r)
#'
#' # 1D PIAB example:
#' num.states <- 4        # Number of PIAB states
#' L.box      <- 50e-6    # Box length
#' m.particle <- 1.78e-37 # Mass of particle(s)
#' TT         <- 77       # Temperature
#' N          <- 2.3*N.A  # Total number of particles
#'
#' # PIAB energies:
#' e <- piab.energy(1:num.states, m = m.particle, L = L.box)
#' e
#'
#--------------------------------------------
piab.energy <- function(n, m, L) {

  # Minor error checking
  if(length(n) != length(L)){
    stop("Length of n and L must be the same!")
  }

  val <- (hb^2 * pi^2)/(2*m) * sum((n^2/L^2))

  return(val)
}


#--------------------------------------------
#' @title Boltzmann factor
#' @description This function computes a Boltzmann factor
#'
#' @param energy an energy, scaled energy or offset energy
#' @param Temp   temperature
#'
#' @details This is a function computes a Boltzmann factor.
#'
#' @return a Boltzmann factor.
#'
#' @examples
#' library(che302r)
#'
#' num.states <- 4        # Number of PIAB states
#' L.box      <- 50e-6    # Box length
#' m.particle <- 1.78e-37 # Mass of particle(s)
#' TT         <- 77       # Temperature
#' N          <- 2.3*N.A  # Total number of particles
#'
#' # PIAB energies:
#' e <- piab.energy(1:num.states, m = m.particle, L = L.box)
#' e
#'
#' # Compute the scaled partition function, Z-tilde:
#' # PIAB energies relative to PIAB state 1:
#' e.tilde <- e-e[1]
#' e.tilde
#'
#' # Compute partition function, Z:
#' # Boltzmann factors for each PIAB state:
#' bfs <- boltz.fac(e, Temp = TT)
#' bfs
#'
#' # Sum of the Boltzmann factors is the partition function
#' Z <- sum(bfs)
#' Z
#'
#' # Boltzmann factors relative to PIAB state 1:
#' bfs.tilde <- boltz.fac(e.tilde, Temp = TT)
#' bfs.tilde
#'
#' # Sum to get Z-tilde:
#' Z.tilde <- sum(bfs.tilde)
#' Z.tilde
#'
#' # Z-tilde repationship to Z
#' Z
#' bfs[1]*Z.tilde # Same as Z??
#'
#--------------------------------------------
boltz.fac <- function(energy, Temp, logQ=F) {

  if(logQ == T) {
    val <- -energy/(kB*Temp)
  } else {
    val <- exp(-energy/(kB*Temp))
  }

  return(val)
}
