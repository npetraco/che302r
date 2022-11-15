#takes wavelength in nm and returns an rgb-alpha value
#From: http://www.physics.sfasu.edu/astro/color/spectra.html converted to R
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

  #   #intensty is lower at the edges of the visible spectrum.
  if (wl > 780 || wl < 380) {
    alpha <- 0
  } else if (wl > 700) {
    alpha <- (780 - wl) / (780 - 700)
  } else if (wl < 420) {
    alpha <- (wl - 380) / (420 - 380)
  } else {
    alpha <- 1
  }

  #colorSpace <- paste("rgba(", (R * 100), "%,", (G * 100), "%,", (B * 100), "%, ", alpha, ")", R, G, B, alpha)
  colorSpace <- rgb(R, G, B, alpha)
  #
  #   #colorSpace is an array with 5 elements.
  #   #The first element is the complete code as a string.
  #   #Use colorSpace[0] as is to display the desired color.
  #   #use the last four elements alone or together to access each of the individual r, g, b and a channels.
  #
  return(colorSpace)

}

#A rough visible spectrum vector:
# wl<-seq(390,730,5)
# rep("purple",18) #390-475
# rep("blue",3) #480-490
# rep("green",13) #495-555
# rep("yellow",4) #560-575
# rep("orange",4) #580-595
rough.spec.vec <- c(rep("purple",8),rep("blue",13),rep("green",13),rep("yellow",4),rep("orange",4),rep("red",27))


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
vec2hex <- function(a.vec){paste(c("#",as.character(as.hexmode(a.vec))), collapse = "")}


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
