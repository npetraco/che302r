#--------------------------------------------
# Quick reference for normalized xy-chromaticity CIE 1931 coordinates in the colorspace library.
# NOTE: needs library(colorspace)
# NOTE: XXXXXX in data
#--------------------------------------------
x31 <- ciexyz31[,2]/(ciexyz31[,2]+ciexyz31[,3]+ciexyz31[,4])
y31 <- ciexyz31[,3]/(ciexyz31[,2]+ciexyz31[,3]+ciexyz31[,4])


#--------------------------------------------
# Spline together CIE 1931 cmf data together so we can use them as function if needed
#--------------------------------------------
xbar.sf <- splinefun(ciexyz31$wlnm, ciexyz31$xbar)
ybar.sf <- splinefun(ciexyz31$wlnm, ciexyz31$ybar)
zbar.sf <- splinefun(ciexyz31$wlnm, ciexyz31$zbar)



#----------------------------------------------------------------
#' @title Generate spectral locus data
#' @description XX
#'
#' @param xch XX
#' @param ych XX
#' @param lambda XX
#'
#' @details XX
#'
#' @return The function will XX
#'
#----------------------------------------------------------------
make.spectral.locus <- function(delta.lambda=1, num.lambda=NULL, lambda.start=380, lambda.stop=780){

  if(is.null(num.lambda)){
    lambda.loc <- seq(from=lambda.start, to=lambda.stop, by=delta.lambda)
  } else {
    lambda.loc <- seq(from=lambda.start, to=lambda.stop, length.out=num.lambda)
  }


  x.spectral.locus <- xbar.sf(lambda.loc)/(xbar.sf(lambda.loc) + ybar.sf(lambda.loc) + zbar.sf(lambda.loc))
  y.spectral.locus <- ybar.sf(lambda.loc)/(xbar.sf(lambda.loc) + ybar.sf(lambda.loc) + zbar.sf(lambda.loc))
  #plot(x.sl, y.sl)

  spectral.locus.mat <- cbind(lambda.loc, x.spectral.locus, y.spectral.locus)
  colnames(spectral.locus.mat) <- c("lambda", "x", "y")

  return(spectral.locus.mat)

}


#----------------------------------------------------------------
#' @title Generate alychne data or line
#' @description XX
#'
#' @param lambda.left XX
#' @param lambda.right XX
#' @param type XX
#' @param num.pts XX
#'
#' @details XX
#'
#' @return The function will XX
#'
#----------------------------------------------------------------
make.alychne <- function(lambda.left=380, lambda.right=780, type="data", num.pts=1000){

  # Left most requested chromaticity point on the alychne
  x.ch.left <- xbar.sf(lambda.left)/(xbar.sf(lambda.left) + ybar.sf(lambda.left) + zbar.sf(lambda.left))
  y.ch.left <- ybar.sf(lambda.left)/(xbar.sf(lambda.left) + ybar.sf(lambda.left) + zbar.sf(lambda.left))

  # Right most requested chromaticity point on the alychne
  x.ch.right <- xbar.sf(lambda.right)/(xbar.sf(lambda.right) + ybar.sf(lambda.right) + zbar.sf(lambda.right))
  y.ch.right <- ybar.sf(lambda.right)/(xbar.sf(lambda.right) + ybar.sf(lambda.right) + zbar.sf(lambda.right))

  m.alychne <- (y.ch.right - y.ch.left)/(x.ch.right-x.ch.left) # alychne slope
  b.alychne <- (y.ch.left-m.alychne*x.ch.left)                         # alychne intercept

  if(type=="data") {
    x.alychne    <- seq(from=x.ch.left, to=x.ch.right, length.out=num.pts)
    y.alychne    <- m.alychne*x.alychne + b.alychne
    alychne.info <- cbind(x.alychne, y.alychne)
    colnames(alychne.info) <- c("x","y")

  } else if(type=="line") {
    alychne.info <- c(m.alychne, b.alychne)
    names(alychne.info) <- c("m", "b")
  } else {
    stop("type must be data or line")
  }

  return(alychne.info)

}


#--------------------------------------------
#' @title wavelengthToColor.
#' @description This is a function takes wavelength in nm and returns an rgb-alpha value in hex.
#'
#' @param wavelength a wavelength in nm
#'
#' @details This is a function takes wavelength in nm and returns an rgb-alpha value in hex.
#' From: http://www.physics.sfasu.edu/astro/color/spectra.html converted to R. I think this
#' guy may have written the original code however: https://earlglynn.github.io/
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

  # intensity is lower at the edges of the visible spectrum.
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


#--------------------------------------------
#' @title Compute "purity" of a point in xy-chromaticity space
#' @description Compute "purity" of a point in xy-chromaticity space
#'
#' @param xc x-chromaticity
#' @param yc y-chromaticity
#' @param xw x-white point. Default is D65 xw=0.31
#' @param yw y-white point. Default is D65 yw=0.33
#' @param percentQ Compute as a percent?
#' @param printQ   Print more info?
#'
#' @details Compute "purity" of a point in xy-chromaticity space  (i.e. "closeness)
#' to the locus convex hull wrt/ a chosen whit point). Assumes as default D65 white point.
#'
#' @return A chromaticity's purity with respect to the white point
#'
#--------------------------------------------
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


#--------------------------------------------
#' @title Compute hue
#' @description Compute hue
#'
#' @param xc x-chromaticity
#' @param yc y-chromaticity
#' @param xw x-white point. Default is D65 xw=0.31
#' @param yw y-white point. Default is D65 yw=0.33
#'
#' @details Compute the approximate closest point on the convex hull of the gamut to
#' chromaticity coords xc, yc wrt white point xw, yw. D65 white point is the
#' default.
#'
#' @return A chromaticity's purity with respect to the white point
#'
#--------------------------------------------
compute_hue <- function(xc, yc, xw=0.31, yw=0.33){

  # defs:
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

  hue.info        <- c(lambda.hue, xh, yh)
  names(hue.info) <- c("lambda", "x", "y")

  return(hue.info)
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

  XYZ.vec.loc <- as.numeric(XYZ.vec)

  T.XYZ2RGB <- rbind(
    c( 3.2406, -1.5372, -0.4986),
    c(-0.9689,  1.8758,  0.0415),
    c( 0.0557, -0.2040,  1.0570)
  )

  RGB.vec <- as.vector(T.XYZ2RGB %*% XYZ.vec.loc)
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

  XYZ.vec.loc <- as.numeric(XYZ.vec)

  T.XYZ2RGB <- rbind(
    c( 2.0413690, -0.5649464, -0.3446944),
    c(-0.9692660,  1.8760108,  0.0415560),
    c( 0.0134474, -0.1183897,  1.0154096)
  )

  RGB.vec <- as.vector(T.XYZ2RGB %*% XYZ.vec.loc)

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

  XYZ.vec.loc <- as.numeric(XYZ.vec)

  T.XYZ2RGB <- rbind(
    c( 2.49349691, -0.93138362, -0.40271078),
    c(-0.82948897,  1.76266406,  0.02362469),
    c( 0.03584583, -0.07617239,  0.95688452)
  )

  RGB.vec <- as.vector(T.XYZ2RGB %*% XYZ.vec.loc)

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


# XYZ2xyY code from colorscience package. I have it in this library because the version
# in colorscience seemed to have a slight bug for certain cases. See code below for the modification made here.
XYZ2xyY.cs <- function (XYZmatrix, illuminant = "D65", observer = 2, RefWhite = get("XYZperfectreflectingdiffuser", envir = environment())) {

  if (!is.matrix(XYZmatrix)) {
    XYZmatrix <- matrix(XYZmatrix, ncol = 3, byrow = TRUE)
  }

  Den <- rowSums(XYZmatrix)
  DenG0 <- which(Den > 0)

  xyYmatrix <- XYZmatrix
  #xyYmatrix[DenG0, 1:2] <- XYZmatrix[DenG0, 1:2]/Den       # BUG?: Original line in colorscience. If some rows in XYZmatrix[DenG0, 1:2] are thrown out because some denominators (Den) are 0, then XYZmatrix[DenG0, 1:2] and Den will have a different number of rows.
  xyYmatrix[DenG0, 1:2] <- XYZmatrix[DenG0, 1:2]/Den[DenG0] # FIX
  xyYmatrix[DenG0, 3]   <- XYZmatrix[DenG0, 2]*100          # Another change; but just cosmetic: Use percent scale for Ys

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
xyY2XYZ2 <- function(xyY.vec){

  x.loc <- xyY.vec[1]
  y.loc <- xyY.vec[2]
  Y.loc <- xyY.vec[3]

  X.loc <- x.loc*Y.loc/y.loc
  Z.loc <- (1-x.loc-y.loc)*Y.loc/y.loc

  return(c(X.loc, Y.loc, Z.loc))
}


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
#----------------------------------------------------------------
plot_spectral_locus <- function(Y.level, white.point=c(0.31,0.33), type="sRGB"){

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

  # hue of each chromaticity on the spectral locus:
  col.nms <- sapply(1:length(cccie31[,1]), function(x){color_partition(cccie31[x,1])})

  #col.info.mat <- data.frame(ciexyz31[,1:3], as.character(col.vec))
  #colnames(col.info.mat) <- c("wlnm","xbar","ybar","col")
  #print(col.info.mat)

  col.info.mat <- data.frame(cccie31[,1], col.nms, col.vec, x31, y31)
  colnames(col.info.mat) <- c("lambda", "hue.partition", "color.code", "x", "y")
  #print(col.info.mat)


  # Identify the wavelengths and chromaticities where the main spectral colors are partitioned:
  # cf. color_partition() for our (rather arbitrary but why not) partition
  col.chunks      <- c("red", "orange", "yellow", "green", "blue", "violet") # ROYGBV
  col.partion.mat <- NULL
  for(i in 1:length(col.chunks)){
    col.hue          <- col.chunks[i]
    col.hue.idxs     <- which(col.info.mat[,2] == col.hue)
    sub.col.info.mat <- col.info.mat[col.hue.idxs,] # Pull out the wl span covering the hue
    col.bottom       <- sub.col.info.mat[1,]        # Bottom of the wl span (i.e. higher energy side)
    #col.top          <- sub.col.info.mat[nrow(sub.col.info.mat),]
    col.partion.mat  <- rbind(col.partion.mat, col.bottom)
  }
  #print(col.partion.mat)

  # For line of purples:
  purples.left       <- col.info.mat[1,]
  purples.right      <- col.info.mat[nrow(col.info.mat),]
  purples.left[2:3]  <- c("purple","#6C3082")  # Call "#6C3082"  purple
  purples.right[2:3] <- c("purple","#6C3082")
  col.partion.mat    <- rbind(col.partion.mat, purples.left, purples.right)
  print(col.partion.mat)

  plot(x31, y31, col=col.vec, pch=16, xlab="x", ylab="y", main="The Great Gamut", xlim=c(0,0.8))
  # Working our way around counter-clockwise:
  lines(c(white.point[1],col.partion.mat[8,4]), c(white.point[2],col.partion.mat[8,5]), col=col.partion.mat[8,2], lwd=3) # designated right purple bound
  lines(c(white.point[1],col.partion.mat[1,4]), c(white.point[2],col.partion.mat[1,5]), col=col.partion.mat[1,2], lwd=3) # designated red bound
  lines(c(white.point[1],col.partion.mat[2,4]), c(white.point[2],col.partion.mat[2,5]), col=col.partion.mat[2,2], lwd=3) # designated orange bound
  lines(c(white.point[1],col.partion.mat[3,4]), c(white.point[2],col.partion.mat[3,5]), col=col.partion.mat[3,2], lwd=3) # designated yellow bound
  lines(c(white.point[1],col.partion.mat[4,4]), c(white.point[2],col.partion.mat[4,5]), col=col.partion.mat[4,2], lwd=3) # designated green bound
  lines(c(white.point[1],col.partion.mat[5,4]), c(white.point[2],col.partion.mat[5,5]), col=col.partion.mat[5,2], lwd=3) # designated blue bound
  lines(c(white.point[1],col.partion.mat[6,4]), c(white.point[2],col.partion.mat[6,5]), col=col.partion.mat[6,2], lwd=3) # designated violet bound
  lines(c(white.point[1],col.partion.mat[7,4]), c(white.point[2],col.partion.mat[7,5]), col=col.partion.mat[7,2], lwd=3) # designated left purple bound
  lines(c(col.partion.mat[7,4],col.partion.mat[8,4]), c(col.partion.mat[7,5],col.partion.mat[8,5]), col=col.partion.mat[7,2], lwd=3) # line of purples

  #return(col.info.dat)

}



#----------------------------------------------------------------
#' @title xy-chromaticity and lambda lookup table
#' @description xy-chromaticity and lambda lookup table
#'
#' @param xch XX
#' @param ych XX
#' @param lambda XX
#'
#' @details xy-chromaticity and lambda lookup table so we can go between xy
#' normalized chromaticities and corresponding lambdas. Have to do this as a
#' lookup table because the gamut is not a function and doesn't have an easy
#' to define inverse.
#'
#' @return The function will XX
#'
#----------------------------------------------------------------
xyl.lut <- function(xch, ych, lambda=NULL, delta.lambda=0.1, tol=0.005, printQ=F){

  # If no lambda is stipulated, check anf see if input xy are on (or close enough too) the gamut convex hull.
  # If they are return the closest xy on the hull to the input ones along with a lambda if they're spectral chromaticities.
  if(is.null(lambda)) {
    spectral.locus <- make.spectral.locus(delta.lambda, lambda.start = 380, lambda.stop = 830) # SHOULD THIS JUST BE A STORED DATA TABLE INSTEAD, SO WE DON'T HAVE TO REGENEATE IT EVERY TIME WE LOOK UP A LAMBDA OR XY??

    dists.from.locus   <- sapply(1:nrow(spectral.locus), function(xx){ed(c(xch, ych), c(spectral.locus[xx,2], spectral.locus[xx,3]))})

    min.dist <- min(dists.from.locus)

    chull.info <- NULL

    # First check to see if xch, ych is closer to the spectral locus or alychne.
    # This is a little sloppy and can be done better.......
    if(min.dist < tol){

      min.idx    <- which(dists.from.locus == min.dist)
      xslocus    <- spectral.locus[min.idx,2]
      yslocus    <- spectral.locus[min.idx,3]
      lamslocus  <- spectral.locus[min.idx,1]
      chull.info <- c(xslocus, yslocus, lamslocus)
      names(chull.info) <- c("x","y","lambda")

      if(printQ==T){
        print("On spectral locus")
        print(paste0("Input chromaticity:            x=", xch, ", y=", ych))
        print(paste0("Closest spectral chromaticity: x=", round(xslocus,3), ", y=", round(yslocus,3), ", lambda=", lamslocus))
        print(paste0("Distance:                       =", min.dist))
      }

    } else { # If not on spectral locus, check to see if xch, ych is close to the alychne using the same distance tolerance:

      alychne.mat        <- make.alychne(lambda.left = 380, lambda.right = 830, type = "data", num.pts = 1000) # SHOULD THIS JUST BE A STORED DATA TABLE INSTEAD, SO WE DON'T HAVE TO REGENEATE IT EVERY TIME WE LOOK UP A LAMBDA OR XY??
      dists.from.alychne <- sapply(1:nrow(alychne.mat), function(xx){ed(c(xch, ych), c(alychne.mat[xx,1], alychne.mat[xx,2]))})
      min.dist.aly       <- min(dists.from.alychne)

      if(min.dist.aly < tol){

        min.aly.idx  <- which(dists.from.alychne == min.dist.aly)
        xaly         <- alychne.mat[min.aly.idx,1]
        yaly         <- alychne.mat[min.aly.idx,2]

        chull.info <- c(xaly, yaly, NA)
        names(chull.info) <- c("x","y", "lambda")

        if(printQ==T){
          print("On alychne")
          print(paste0("Input chromaticity:            x=", xch, ", y=", ych))
          print(paste0("Closest alychne chromaticity:  x=", round(xaly,3), ", y=", round(yaly,3)))
          print(paste0("Distance:                       =", min.dist.aly))
        }

      } else {

        min.aly.idx  <- which(dists.from.alychne == min.dist.aly)
        print("Chromaticity not on spectral locus or alychne!")
        print(paste0(xch, " ", ych))
        print(paste0(round(alychne.mat[min.aly.idx,1], 3), " ", round(alychne.mat[min.aly.idx,2], 3) ))
      }

    }

  } else { # If a lambda is stipulated compute corresponding spectral xy and ignore input xy if they're there:

    xslocus           <- xbar.sf(lambda)/(xbar.sf(lambda) + ybar.sf(lambda) + zbar.sf(lambda))
    yslocus           <- ybar.sf(lambda)/(xbar.sf(lambda) + ybar.sf(lambda) + zbar.sf(lambda))
    chull.info        <- c(xslocus, yslocus, lambda)
    names(chull.info) <- c("x","y","lambda")

  }

  return(chull.info)

}
