#--------------------------------------------
#' @title Constants
#'
#' @description Constants available
#'
#' @details  These physical constants and handy conversions are available in the che302r library.
#' Reference for these and more cf. https://www.nist.gov/pml/fundamental-physical-constants
#'
#' @usage Constants:
#' h      Planck's const, Js
#' cl     Speed of light, m/s
#' kB     Boltzmann's const, J/K
#' hb     Reduced Planck's const, Js
#' N.A    Avogadro's number, particles
#' me     mass of electron, kg
#' mP     mass of proton, kg
#' mN     mass of neutron, kg
#' ec     charge of electron, Coulombs
#' RH     Rydberg constant for hydrogen atom, wavenumbers cm^-1
#' RHJ    Rydberg constant for hydrogen atom, Joules
#' Rinf   Rydberg constant clamped nucleus, wavenumbers cm^-1
#' RinfJ  Rydberg constant clamped nucleus, Joules
#' alphaf Fine structure constant, dimensionless
#' eps0   Vacuum permittivity, Farads/m
#' a0     Bohr radius, m
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
#' @export
#--------------------------------------------
constants <- function(){
  print("See: ?constants")
}

h      <- 6.62607015e-34           # Planck's const NIST, Js
cl     <- 299792458                # Speed of light NIST, m/s
kB     <- 1.380649e-23             # Boltzmann's const NIST, J/K
kb     <- kB                       # Also Boltzmann's const, J/K
hb     <- h/(2*pi)                 # Reduced Planck's const, Js
hbar   <- hb                       # Also reduced Planck's const, Js
h.bar  <- hb                       # Also^2 reduced Planck's const, Js
N.A    <- 6.02214076e23            # Avogadro's number NIST, particles
na     <- N.A                      # Also Avogadro's number, particles; against my better judgement......
me     <- 9.1093837139e-31         # Mass of electron, NIST, kg
mP     <- 1.67262192595e-27        # Mass of proton NIST, kg
mN     <- 1.67492750056e-27        # Mass of neutron NIST, kg
muH    <- me*mP/(me+mP)            # Reduced mass of Hydrogen atom
ec     <- 1.602176634e-19          # Charge of electron from NIST, Coulombs
alphaf <- 0.0072973525643          # Fine structure constant from NIST
eps0   <- 8.8541878188e-12         # Vacuum permittivity NIST, Farads/m
eps0c  <- ec^2/(2*alphaf*h*cl)     # Vacuum permittivity from constants, Farads/m
a0     <- 5.29177210544e-11        # Bohr radius from NIST, m
a0c    <- 4*pi*eps0*hb^2/(me*ec^2) # Bohr radius from constants
Rinf       <- 10973731.568157               # Rinfinity Rydberg constant from NIST, wavenumbers m^-1
RinfJ.nist <- 2.1798723611030e-18           # Rinfinity Rydberg constant in Joules from NIST
Rinfc      <- alphaf/(4*pi*a0)              # Rinfinity Rydberg constant from constants, wavenumbers m^-1
RinfJ      <- h * cl/(1/(Rinf))             # Rinfinity Rydberg constant NIST value for Rinfinity, Joules
RinfJc     <- me*ec^4/(32*pi^2*eps0^2*hb^2) # Rinfinity Rydberg constant from other constants, Joules
RH         <- muH/me * Rinf/100             # "Corrected" Rydberg constant for Hydrogen, wavenumbers cm^-1
RHJ        <- h * cl/(1/(RH*100))           # "Corrected" Rydberg constant for Hydrogen converted from RH, Joules

# Handy Conversions
eV2J         <- 1.602176634e-19                         # electron-volt to Joules, energy gained by 1 e- accelerated though 1V, i.e. 1ec*1V. Same magnitude as ec bec 1V = 1J/C
J2eV         <- eV2J^(-1)                               # Joules to eV
hartree2J    <- 4.3597447222060e-18                     # hartrees to Joules from NIST
amu2kg       <- 1.0/(N.A*1000)                          # amu to kg = 1g/mol * 1mol/N.A * 1kg/1000g
kg2amu       <- N.A*1000                                # kg to amu
bohr2m       <- a0                                      # bohr to meters
auf2Hz       <- sqrt(hartree2J * 1/(amu2kg * bohr2m^2)) # atomic units frequency to Hz, Unit= s^-1
lambdaFw2icm <- auf2Hz/(2*pi*cl*100)                    # atomic units mass weighted Fw eigenvalues to nu-tilde = omega/(2*pi*c), Unit = cm^-1
                                                        # lambdaF2icm ---> 5140.485 # Conversion factor from sqrt(hartree/amu-bohr^2) to cm^-1

#--------------------------------------------
#' @title "Corrected" Rydberg constant for any hydrogen like atom
#' @description "Corrected" Rydberg constant for any hydrogen like atom
#'
#' @param at.mass   atomic nucleus' mass
#' @param mass.unit "amu" or "kg"
#' @param en.unit   Rydberg constant unit desired: "cm-1" (inverse centimeters), "m-1" (inverse meters) or J ("Joules)
#'
#' @details The Rinf (or RinfJ) Rydberg constant assumes the center of mass for an nucleus-electron system is
#' directly at the nucleus (i.e. nucleus is clamped). To correct for this multiply by mu = me*mnuc/(me+mnuc)
#' and divide by me.
#'
#' @return The corrected hydrogen-like atom Rydberg constant
#'
#' @examples
#' RydM(mP, mass.unit = "kg", en.unit = "cm-1")
#' RH
#'
#' RydM(1, mass.unit = "amu", en.unit = "cm-1")
#' RH # a tad different because 1 amu != mP or mN
#'
#' RydM(2*mP+2*mN, mass.unit = "kg", en.unit = "cm-1") # Rydberg for Helium-4 +1 ion
#--------------------------------------------
RydM <- function(at.mass, mass.unit="amu", en.unit="cm-1"){

  if(mass.unit == "amu"){
    avg.nuc.mass.kg <- at.mass * amu2kg # For amu masses assume average isotope mass in kg and convert using amu2kg
    mu.at <- avg.nuc.mass.kg*me/(avg.nuc.mass.kg + me)
  } else if(mass.unit =="kg"){
    mu.at <- at.mass*me/(at.mass + me)
  } else {
    stop("mass.unit must be amu or kg!")
  }

  if(en.unit == "m-1"){ # Unit here is m^-1
    val <- mu.at/me*Rinf
  } else if(en.unit == "cm-1") { # Unit here is cm^-1
    val <- mu.at/me*(Rinf/100)
  } else if(en.unit == "J") { # Unit here is J
    val <- mu.at/me*RinfJ
  } else {
    stop("en.unit must be cm-1, m-1 or J")
  }

  return(val)

}
