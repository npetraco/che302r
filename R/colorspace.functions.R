#--------------------------------------------
#' @title wavelengthToColor.
#' @description This is a function takes wavelength in nm and returns an rgb-alpha value in hex.
#'
#' @param wavelength a wavelength in nm
#'
#' @details This is a function takes wavelength in nm and returns an rgb-alpha value in hex.
#' From: http://www.physics.sfasu.edu/astro/color/spectra.html converted to R.
#'
#' @return A hexcode for an RGB color
#'
#' @examples
#' wavelengthToColor(536)
#--------------------------------------------
wavelengthToColor <- function(wavelength) {
  wl <- wavelength
  gamma <- 1

  if (wl >= 380 && wl < 440) {
    R <- -1 * (wl - 440) / (440 - 380)
    G <- 0
    B <- 1
  } else if (wl >= 440 && wl < 490) {
    R <- 0
    G <- (wl - 440) / (490 - 440)
    B <- 1
  } else if (wl >= 490 && wl < 510) {
    R <- 0
    G <- 1
    B <- -1 * (wl - 510) / (510 - 490)
  } else if (wl >= 510 && wl < 580) {
    R <- (wl - 510) / (580 - 510)
    G <- 1
    B <- 0
  } else if (wl >= 580 && wl < 645) {
    R <- 1
    G <- -1 * (wl - 645) / (645 - 580)
    B <- 0.0
  } else if (wl >= 645 && wl <= 780) {
    R <- 1
    G <- 0
    B <- 0
  } else {
    R <- 0
    G <- 0
    B <- 0
  }

  # intensty is lower at the edges of the visible spectrum.
  if (wl > 780 || wl < 380) {
    alpha <- 0
  } else if (wl > 700) {
    alpha <- (780 - wl) / (780 - 700)
  } else if (wl < 420) {
    alpha <- (wl - 380) / (420 - 380)
  } else {
    alpha <- 1
  }

  # colorSpace is an array with 5 elements.
  # The first element is the complete code as a string.
  # Use colorSpace[0] as is to display the desired color.
  # Use the last four elements alone or together to access
  # each of the individual r, g, b and a channels.
  #
  #colorSpace <- paste("rgba(", (R * 100), "%,", (G * 100), "%,", (B * 100), "%, ", alpha, ")", R, G, B, alpha)

  # Do it this way instead of the above:
  colorSpace <- rgb(R, G, B, alpha)

  return(colorSpace)
}


#--------------------------------------------
#' @title color_partition
#' @description This is a function takes wavelength in nm and returns simple color name
#'
#' @param wavelength a wavelength in nm
#'
#' @details This is a function takes wavelength in nm and returns simple color name.
#' It works the same way as wavelengthToColor
#'
#' @return A hexcode for an RGB color
#'
#' @examples
#' wavelengthToColor(536)
#--------------------------------------------
color_partition <- function(wavelength) {
  wl <- wavelength

  if (wl >= 380 && wl < 440) {
    color.name <- "violet"
  } else if (wl >= 440 && wl < 490) {
    color.name <- "blue"
  } else if (wl >= 490 && wl < 510) {
    color.name <- "green"
  } else if (wl >= 510 && wl < 580) {
    color.name <- "yellow"
  } else if (wl >= 580 && wl < 645) {
    color.name <- "orange"
  } else if (wl >= 645 && wl <= 780) {
    color.name <- "red"
  } else {
    #color.name <- NA
    color.name <- "#00000000"
  }

  return(color.name)

}


#A rough visible spectrum vector:
# wl<-seq(390,730,5)
# rep("purple",18) #390-475
# rep("blue",3) #480-490
# rep("green",13) #495-555
# rep("yellow",4) #560-575
# rep("orange",4) #580-595
rough.spec.vec <- c(rep("purple",8),rep("blue",13),rep("green",13),rep("yellow",4),rep("orange",4),rep("red",27))

#normalized xy chromaticity coordinates. NOTE: needs library(colorspace)
x31 <- ciexyz31[,2]/(ciexyz31[,2]+ciexyz31[,3]+ciexyz31[,4])
y31 <- ciexyz31[,3]/(ciexyz31[,2]+ciexyz31[,3]+ciexyz31[,4])


# Compute "purity" of a point in xy-chromaticity space (i.e. it's "closeness"
# to the locus convex hull wrt/ a chosen whit point). Assumes as default D65 white point.
color_purity <- function(xc, yc, xw=0.31, yw=0.33, percentQ=T, printQ=F){

  # Compute the nearest point on the gamut convex hull to xc, yc wrt/ the white point:
  hue.info.vec <- compute_hue(xc=xc, yc=yc, xw=xw, yw=yw)
  hue.lambda   <- hue.info.vec[1]
  xhull        <- hue.info.vec[2]
  yhull        <- hue.info.vec[3]

  val <- sqrt(((xc-xw)^2 +(yc-yw)^2)/((xhull-xw)^2 +(yhull-yw)^2))
  if(percentQ==T){
    val <- val*100
  }

  if(printQ==T){
    print(paste0("Comparing to lambda: ", hue.lambda))
    print(paste0("Lambda is a:         ", color_partition(hue.lambda)))
    if(percentQ==T){
      print(paste0("Purity (saturation): ", round(val, 1), "%"))
    } else {
      print(paste0("Purity (saturation): ", round(val, 1)))
    }
  }

  return(val)

}


# Compute the approximate closest point on the convex hull of the gamut to
# chromaticity coords xc, yc wrt white point xw, yw. D65 white point is the
# default.
compute_hue <- function(xc, yc, xw=0.31, yw=0.33){

  # defs.
  # chromaticity vec: vector from white point to (xc,yc) chromaticity value
  # spectral vec: vector from white point to a point on the spectral locus (edge of the gamut)

  cxy <- c(xc,yc)
  wp  <- c(xw,yw)
  rxy <- cxy - wp

  angs <- rep(0, nrow(ciexyz31))
  prjs <- rep(0, nrow(ciexyz31))
  for(i in 1:nrow(ciexyz31)){
    rs      <- c(x31[i],y31[i]) - wp
    prjs[i] <- rxy%*%rs/sqrt(sum(rs^2))                # projection of chromaticity vec onto spectral vec
    angs[i] <- acos(prjs[i]/sqrt(sum(rxy^2))) * 180/pi # compute angle between chromaticity vec and spectral vec
  }

  #print(cbind(x31,y31,prjs,angs))
  #return(cbind(x31,y31,prjs,angs))

  mind.idx   <- which(angs==min(angs)) # Choose the point that makes the smallest angle with the chromaticity vector. NOTE: the chosen point lies on the spectral vector that makes the smallest angle with the chromaticity vector
  lambda.hue <- ciexyz31[mind.idx,1]
  xh         <- x31[mind.idx]
  yh         <- y31[mind.idx]

  return(c(lambda.hue, xh, yh))
}


# Gamma correction for sRGB
gamma.sRGB <- function(RGB.vec){
  RGB.vec.gc <- RGB.vec
  for(i in 1:3){
    if(RGB.vec[i] <= 0.0031308){
      RGB.vec.gc[i] <-12.92*RGB.vec[i]
    } else {
      RGB.vec.gc[i] <- 1.055*(RGB.vec[i])^(1/2.4) - 0.055
    }
  }
  return(RGB.vec.gc)
}


# Gamma correction for AdobeRGB
gamma.AdobeRGB <- function(RGB.vec){

  gammaa     <- 1/2.19921875
  RGB.vec.gc <- RGB.vec

  for(i in 1:3){
    if(RGB.vec[i] >= 0){
      RGB.vec.gc[i] <-RGB.vec[i]^gammaa
    } else {
      RGB.vec.gc[i] <- -1 * (-1*RGB.vec[i])^gammaa
    }
  }
  return(RGB.vec.gc)

}


# Gamma correction for P3RGB
gamma.P3RGB <- function(RGB.vec){

  gammaa     <- 1/2.6 # ????
  RGB.vec.gc <- RGB.vec

  for(i in 1:3){
    if(RGB.vec[i] >= 0){
      RGB.vec.gc[i] <-RGB.vec[i]^gammaa
    } else {
      RGB.vec.gc[i] <- -1 * (-1*RGB.vec[i])^gammaa
    }
  }
  return(RGB.vec.gc)

}


# Convert a vector to hex components formatted conveniently to use as a color code
vec2hex <- function(a.vec){

  a.vec.hex <- NULL
  for(i in 1:length(a.vec)){
    if(abs(a.vec[i]) > 15) {
      hex1i     <- as.character(as.hexmode(a.vec[i]))
      a.vec.hex <- paste(a.vec.hex, hex1i, sep="")
    } else if(abs(a.vec[i]) <= 15) {
      # as.hexmode() just gives four bits back for numbers 0-15. Paste on an extra 0 needed for color conversion
      hex1i <- as.character(as.hexmode(a.vec[i]))
      a.vec.hex <- paste(a.vec.hex, "0", hex1i, sep="")
    } else {
      stop("Vector's components must be numbers!")
    }
  }
  a.vec.hex <- paste0("#",a.vec.hex)

  return(a.vec.hex)
}


#----------------------------------------------------------------
#' @title Euclidean distance
#' @description Euclidean distance in a plane
#'
#' @param vec1 a vector
#' @param vec2 a vector
#'
#' @details Euclidean distance in a plane
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
ed <- function(vec1, vec2){
  return(sqrt(sum((vec1 - vec2)^2)))
}


#----------------------------------------------------------------
#' @title Magnitude of a vector
#' @description Magnitude of a vector
#'
#' @param vec1 a vector
#' @param vec2
#'
#' @details Magnitude of a vector. Handy for working in the gamut
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
magv <- function(vec1){
  return(sqrt(sum(vec1^2)))
}


#----------------------------------------------------------------
#' @title Scalar projection of vector 1 on to vector 2
#' @description Scalar projection of vector 1 on to vector 2
#'
#' @param vec1 vector 1
#' @param vec2 vector 2
#'
#' @details Scalar projection of vector 1 on to vector 2. Usually vector 1 is the
#' shorter of the two vectors, and vector 2 is the longer.
#' Handy for working in the gamut.
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
prjv <- function(vec1,vec2){
  return(vec1%*%vec2/magv(vec2))
}


#----------------------------------------------------------------
#' @title Angle between two vectors
#' @description Angle between two vectors
#'
#' @param vec1 vector 1
#' @param vec2 vector 2
#'
#' @details Angle between two vectors. Handy for working in the gamut
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
angv <- function(vec1,vec2,type="deg"){

  theta <- acos( prjv(vec1,vec2)/magv(vec1) )
  if(type=="deg"){
    theta <-  theta * 180/pi
  }

  return(theta)
}


#----------------------------------------------------------------
#' @title Area for a triangle by 3x3 determinant of vertices
#' @description Area for a triangle by 3x3 determinant of vertices
#'
#' @param vtx1 vertex 1
#' @param vtx2 vertex 2
#' @param vtx3 vertex 3
#'
#' @details Area for a triangle by 3x3 determinant of vertices.
#' Handy for testing if chromaticity point is in purple region.
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
tria <- function(vtx1, vtx2, vtx3){

  tmat <- rbind(
    c(vtx1,1),
    c(vtx2,1),
    c(vtx3,1)
  )

  tarea <- abs(0.5*det(tmat))

  return(tarea)
}


#----------------------------------------------------------------
#' @title Test to see if a point in the gamut is in the purple region
#' @description Test to see if a point in the gamut is in the purple region
#'
#' @param cp                  Chromaticity point in the gamut
#' @param white.point         A white point. Default is D65
#' @param left.purple.vertex  A left purple vertex. Default is the one in ciexyz31
#' @param right.purple.vertex A right purple vertex. Default is the one in ciexyz31
#'
#' @details Test to see if a point in the gamut is in the purple region. The purple
#' region is a convex hull.
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
purpleQ <- function(cp, white.point=c(0.31,0.33), left.purple.vertex=c(0.175560232,0.005293837), right.purple.vertex=c(0.73469,0.26531), tol=1e-5){

  # Calculate area of the whole purple region triangle: lp-wp-rp
  AA <- tria(left.purple.vertex, white.point, right.purple.vertex)

  # Area of triangle: cp-wp-lp
  BB <- tria(cp, white.point, left.purple.vertex)

  # Area of triangle: cp-lp-rp
  CC <- tria(cp, left.purple.vertex, right.purple.vertex)

  # Area of triangle: cp-wp-rp
  DD <- tria(cp, white.point, right.purple.vertex)

  # If the area of the whole triangle is the sum of the areas of the partitions (to tol),
  # then cp is in the purple region.
  if (abs(BB+CC+DD - AA) < tol) {
    inpurplesQ <- TRUE
  } else {
    inpurplesQ <- FALSE
  }

  return(inpurplesQ)
}


# Hand convert to sRGB:
XYZ2sRGB <- function(XYZ.vec){

  T.XYZ2RGB <- rbind(
    c( 3.2406, -1.5372, -0.4986),
    c(-0.9689,  1.8758,  0.0415),
    c( 0.0557, -0.2040,  1.0570)
  )

  RGB.vec <- as.vector(T.XYZ2RGB %*% XYZ.vec)
  #print("Step1")
  #print(RGB.vec)

  # gamma correction:
  RGB.vec <- gamma.sRGB(RGB.vec)
  #print("Step2")
  #print(RGB.vec)

  # Re-scale from 0 to 255
  RGB.vec <- round(RGB.vec*255)
  #print("Step3")
  #print(RGB.vec)

  #Clip any negatives to 0
  for(i in 1:3){
    if(RGB.vec[i] < 0) {
      RGB.vec[i] <- 0
    }
  }
  #print("Step4")
  #print(RGB.vec)

  #Clip any +255 to 255
  for(i in 1:3){
    if(RGB.vec[i] > 255) {
      RGB.vec[i] <- 255
    }
  }
  #print("Step5")
  #print(RGB.vec)

  # Convert to a hex code:
  RGB.hex <- vec2hex(RGB.vec)
  #print("Step6")
  #print(RGB.hex)

  # Check and see if RGB all ended up as 0......
  if(sum(RGB.vec) == 0) {
    RGB.hex <- "#000000"
    warning("All three channels = 0!")
  }

  RGB.info         <- list(RGB.vec, RGB.hex)
  names(RGB.info)  <- c("RGB","hex")

  return(RGB.info)

}


# Hand convert to AdobeRGB:
XYZ2Adobe <- function(XYZ.vec){

  T.XYZ2RGB <- rbind(
    c( 2.0413690, -0.5649464, -0.3446944),
    c(-0.9692660,  1.8760108,  0.0415560),
    c( 0.0134474, -0.1183897,  1.0154096)
  )

  RGB.vec <- as.vector(T.XYZ2RGB %*% XYZ.vec)

  # gamma correction:
  RGB.vec <- gamma.AdobeRGB(RGB.vec)

  # Re-scale from 0 to 255
  RGB.vec <- round(RGB.vec*255)
  #print("Step3")
  #print(RGB.vec)

  #Clip any negatives to 0
  for(i in 1:3){
    if(RGB.vec[i] < 0) {
      RGB.vec[i] <- 0
    }
  }
  #print("Step4")
  #print(RGB.vec)

  #Clip any +255 to 255
  for(i in 1:3){
    if(RGB.vec[i] > 255) {
      RGB.vec[i] <- 255
    }
  }
  #print("Step5")
  #print(RGB.vec)

  # Convert to a hex code:
  RGB.hex <- vec2hex(RGB.vec)
  #print("Step6")
  #print(RGB.hex)

  # Check and see if RGB all ended up as 0......
  if(sum(RGB.vec) == 0) {
    RGB.hex <- "#000000"
    warning("All three channels = 0!")
  }

  RGB.info         <- list(RGB.vec, RGB.hex)
  names(RGB.info)  <- c("RGB","hex")

  return(RGB.info)

}


# Hand convert to DCI-P3 Display:
XYZ2P3 <- function(XYZ.vec){

  T.XYZ2RGB <- rbind(
    c( 2.49349691, -0.93138362, -0.40271078),
    c(-0.82948897,  1.76266406,  0.02362469),
    c( 0.03584583, -0.07617239,  0.95688452)
  )

  RGB.vec <- as.vector(T.XYZ2RGB %*% XYZ.vec)

  # gamma correction:
  RGB.vec <- gamma.P3RGB(RGB.vec)

  # Re-scale from 0 to 255
  RGB.vec <- round(RGB.vec*255)
  #print("Step3")
  #print(RGB.vec)

  #Clip any negatives to 0
  for(i in 1:3){
    if(RGB.vec[i] < 0) {
      RGB.vec[i] <- 0
    }
  }
  #print("Step4")
  #print(RGB.vec)

  #Clip any +255 to 255
  for(i in 1:3){
    if(RGB.vec[i] > 255) {
      RGB.vec[i] <- 255
    }
  }
  #print("Step5")
  #print(RGB.vec)

  # Convert to a hex code:
  RGB.hex <- vec2hex(RGB.vec)
  #print("Step6")
  #print(RGB.hex)

  # Check and see if RGB all ended up as 0......
  if(sum(RGB.vec) == 0) {
    RGB.hex <- "#000000"
    warning("All three channels = 0!")
  }

  RGB.info         <- list(RGB.vec, RGB.hex)
  names(RGB.info)  <- c("RGB","hex")

  return(RGB.info)

}


#----------------------------------------------------------------
#' @title Plot a color with specified chromaticity and brightness values
#' @description Plot a color with specified chromaticity and brightness values
#'
#' @param xyY.vex an x, y, Y color specification
#' @param type    subgamut to use: Adobe, sRGB or P3
#'
#' @details Plot a color with specified chromaticity and brightness values
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
xyY_swatch <- function(xyY.vec, type="Adobe"){

  x.loc <- xyY.vec[1]
  y.loc <- xyY.vec[2]
  Y.loc <- xyY.vec[3]

  X.loc <- x.loc*Y.loc/y.loc
  Z.loc <- (1-x.loc-y.loc)*Y.loc/y.loc

  if(type=="Adobe"){
    hex.cod <- XYZ2Adobe(c(X.loc, Y.loc, Z.loc))$hex
  } else if(type=="sRGB") {
    hex.cod <- XYZ2sRGB(c(X.loc, Y.loc, Z.loc))$hex
  } else if(type=="P3") {
    hex.cod <- XYZ2P3(c(X.loc, Y.loc, Z.loc))$hex
  } else {
    stop("type must be Adobe, sRGB or P3")
  }

  plot(1, pch=16, cex=16, col=hex.cod)

  return(hex.cod) # return it if we want to use it somewhere else

}


# XYZ2xyY from colorscience package. I have it in this library because the version
# in colorscience seemed to have a slight bug. See code below for the modification made here.
XYZ2xyY.cs <- function (XYZmatrix, illuminant = "D65", observer = 2, RefWhite = get("XYZperfectreflectingdiffuser", envir = environment())) {

  if (!is.matrix(XYZmatrix)) {
    XYZmatrix <- matrix(XYZmatrix, ncol = 3, byrow = TRUE)
  }

  Den <- rowSums(XYZmatrix)
  DenG0 <- which(Den > 0)

  xyYmatrix <- XYZmatrix
  #xyYmatrix[DenG0, 1:2] <- XYZmatrix[DenG0, 1:2]/Den       # BUG?: Original line in colorscience. If some rows are thrown out because some denominators (Den) are 0, then XYZmatrix[DenG0, 1:2] and Den will have a different number of rows.
  xyYmatrix[DenG0, 1:2] <- XYZmatrix[DenG0, 1:2]/Den[DenG0] # FIX
  xyYmatrix[DenG0, 3] <- XYZmatrix[DenG0, 2]*100            # Cosmetic change: Use percent scale for Ys

  R <- RefWhite[which(RefWhite[["Illuminant"]] == illuminant),]
  Rx <- unlist(R[paste("X", observer, sep = "")])
  Ry <- unlist(R[paste("Y", observer, sep = "")])
  Rz <- unlist(R[paste("Z", observer, sep = "")])
  x <- Rx/(Rx + Ry + Rz)
  y <- Ry/(Rx + Ry + Rz)

  xyYmatrix[-DenG0, 1] <- x
  xyYmatrix[-DenG0, 2] <- y

  return(xyYmatrix)
}


# An adequate version exists in colorscience package, so comment out for now.
# Convenience function to convert xyY to XYZ
# xyY2XYZ <- function(xyY.vec){
#
#   x.loc <- xyY.vec[1]
#   y.loc <- xyY.vec[2]
#   Y.loc <- xyY.vec[3]
#
#   X.loc <- x.loc*Y.loc/y.loc
#   Z.loc <- (1-x.loc-y.loc)*Y.loc/y.loc
#
#   return(c(X.loc, Y.loc, Z.loc))
# }


#----------------------------------------------------------------
#' @title Plot the outline of the gamut with and without partitions
#' @description Plot the outline of the gamut with and without partitions
#'
#' @param Y.level     The brightness value (Y) to use
#' @param partitionsQ Whether or not to show partitions between the rough hues
#' @param white.point White point. Default is D65
#' @param type        Color sub-gamut to use, Adobe, sRGB, P3, rough.hue, wl2c
#'
#' @details Plot the outline of the gamut with and without partitions
#'
#' @return The function will XX
#'
#' @export
#----------------------------------------------------------------
plot_spectral_locus <- function(Y.level, white.point=c(0.31,0.33), type="Adobe"){

  if(type=="Adobe") {
    col.vec <- sapply(1:nrow(ciexyz31), function(xx){XYZ2Adobe(xyY2XYZ(c(x31[xx], y31[xx], Y.level)))$hex})
  } else if(type=="sRGB") {
    col.vec <- sapply(1:nrow(ciexyz31), function(xx){XYZ2sRGB(xyY2XYZ(c(x31[xx], y31[xx], Y.level)))$hex})
  } else if(type=="P3") {
    col.vec <- sapply(1:nrow(ciexyz31), function(xx){XYZ2P3(xyY2XYZ(c(x31[xx], y31[xx], Y.level)))$hex})
  } else if(type=="rough.hue") {
    col.vec <- sapply(1:nrow(ciexyz31), function(xx){color_partition(wavelength = ciexyz31[xx,1])})
  } else if(type=="wl2c") {
    col.vec <- sapply(1:nrow(ciexyz31), function(xx){wavelengthToColor(wavelength = ciexyz31[xx,1])})
  } else {
    stop("type must be Adobe, sRGB, P3, rough.hue, wl2c")
  }

  col.info.dat <- data.frame(ciexyz31[,1:3], as.character(col.vec))
  colnames(col.info.dat) <- c("wlnm","xbar","ybar","col")



  return(col.info.dat)

}
