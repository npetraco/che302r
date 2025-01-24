#' @title Constants
#'
#' @description Constants available
#'
#' @details  These physical constants and handy conversions are avilible in the che302r library.
#'
#' @usage Constants:
#' h   Planck's const, Js
#' cl  Speed of light, m/s
#' kB  Boltzmann's const, J/K
#' hb  Reduced Planck's const, Js
#' N.A Avogadro's number, particles
#' me  mass of electron, kg
#' mP  mass of proton, kg
#' mN  mass of neutron, kg
#' ec  charge of electron, Coulombs
#' RH  Rydberg constant, wavenumbers cm^-1
#' RHJ Rydberg constant, Joules
#'
#' Handy Conversions:
#' eV2J         electron-volt to Joules, energy gained by 1 e- accelerated though 1V
#' J2eV         Joules to eV
#' hartree2J    hartrees to Joules
#' amu2kg       amu to kg
#' bohr2m       bohr to meters
#' auf2Hz       atomic units frequency to Hz, Unit= s^-1
#' lambdaFw2icm atomic units mass weighted Fw eigenvalues to nu-tilde = omega/(2*pi*c), Unit = cm^-1
#' lambdaF2icm  Conversion factor from sqrt(hartree/amu-bohr^2) to cm^-1
constants <- function(){
  print("See: ?constants")
}
h     <- 6.62607015e-34        # Planck's const, Js
cl    <- 299792458             # Speed of light, m/s
kB    <- 1.3806503e-23         # Boltzmann's const, J/K
kb    <- 1.3806503e-23         # Also Boltzmann's const, J/K
hb    <- h/(2*pi)              # Reduced Planck's const, Js
hbar  <- h/(2*pi)              # Also reduced Planck's const, Js
h.bar <- h/(2*pi)              # Also^2 reduced Planck's const, Js
N.A   <- 6.02214076e23         # Avogadro's number, particles
na    <- N.A                   # Also Avogadro's number, particles; against my better judgement......
me    <- 9.1093837015e-31      # mass of electron, kg
mP    <- 1.67262192369e-27     # mass of proton, kg
mN    <- 1.67492749804e-27     # mass of neutron, kg
ec    <- 1.60217663e-19        # charge of electron, Coulombs
RH    <- 109625                # Rydberg constant, wavenumbers cm^-1
RHJ   <- h * cl/( 1/(RH*100) ) # Rydberg constant, Joules

# Handy Conversions
eV2J         <- 1.602176634e-19                         # electron-volt to Joules, energy gained by 1 e- accelerated though 1V
J2eV         <- eV2J^(-1)                               # Joules to eV
hartree2J    <- 4.3597482e-18                           # hartrees to Joules
amu2kg       <- 1.6603145e-27                           # amu to kg
bohr2m       <- 5.2918e-11                              # bohr to meters
auf2Hz       <- sqrt(hartree2J * 1/(amu2kg * bohr2m^2)) # atomic units frequency to Hz, Unit= s^-1
lambdaFw2icm <- auf2Hz/(2*pi*cl*100)                    # atomic units mass weighted Fw eigenvalues to nu-tilde = omega/(2*pi*c), Unit = cm^-1
                                                        # lambdaF2icm ---> 5140.485 # Conversion factor from sqrt(hartree/amu-bohr^2) to cm^-1
