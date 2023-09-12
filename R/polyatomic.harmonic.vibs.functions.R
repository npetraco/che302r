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
