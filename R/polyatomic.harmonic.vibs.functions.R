hartree2J    <- 4.3597482e-18                           # hartrees to Joules
amu2kg       <- 1.6603145e-27                           # amu to kg
bohr2m       <- 5.2918e-11                              # bohr to meters
auf2Hz       <- sqrt(hartree2J * 1/(amu2kg * bohr2m^2)) # atomic units frequency to Hz, Unit= s^-1
lambdaFw2icm <- auf2Hz/(2*pi*cl*100)                    # atomic units mass weighted Fw eigenvalues to nu-tilde = omega/(2*pi*c), Unit = cm^-1
                                                        #lambdaF2icm ---> 5140.485 # Conversion factor from sqrt(hartree/amu-bohr^2) to cm^-1

#--------------------------------------------
#' @title XXXX
#' @description XXXX
#'
#' @param nqn XXXX
#'
#' @details XXXX
#'
#' @return XXXXX
#'
#' @references XXXXX
#'
#' @examples
#' library(che302r)
#' XXXX
#--------------------------------------------
atom2color<-function(atom.char) {

  atom.colors <- c(
    "H"="green",
    "C"="black",
    "O"="red",
    "N"="blue"
    )

  atom.col <- atom.colors[atom.char]
  return(atom.col)

}

#--------------------------------------------
#' @title XXXX
#' @description XXXX
#'
#' @param nqn XXXX
#'
#' @details XXXX
#'
#' @return XXXX
#'
#' @references XXXXX
#'
#' @examples
#' library(che302r)
#' XXXX
#--------------------------------------------
atom.colors<-function(atoms.char.vec) {

  col.vec <- sapply(1:length(atoms.char.vec), function(xx){atom2color(atoms.char.vec[xx])})
  return(col.vec)

}


#--------------------------------------------
#' @title XXXX
#' @description XXXX
#'
#' @param nqn XXXX
#'
#' @details XXXX
#'
#' @return XXXX
#'
#' @references XXXXX
#'
#' @examples
#' library(che302r)
#' XXXX
#--------------------------------------------
oscillation.f<-function(xyz0.coords, Qj.mat, t.coord) {

  xyzt <- xyz0.coords + Qj.mat*sin(t.coord)

  return(xyzt)
}
