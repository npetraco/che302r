#Planck Black Body Distribution:
#nt = nu-tilde (m^-1)
#Temp = temperature (K)
B<-function(nt,Temp) {
  val<- (2*h*cl^2*nt^3* 1/(exp((h*cl*nt)/(kB*Temp)) - 1) )
  return(val)
}

#Function to generate an interferogram as a function of mirror position in an interferometer
#Input a wavenumber range in units of m^-1 and the corresponding intensity distribution
#wns = wave number range
#Bp = corresponding intensity distribution
#sft = mirror position
#NOTE: THIS IS SLOW AND IS TEMEPERMENTAL AROUND 0 m^-1
interferogram.func<-function(wns,Bp,sft) {
  integrand.func<-splinefun(wns,Bp*cos(2*pi*wns*sft))
  val<-integrate(integrand.func,lower=min(wns),upper=max(wns),stop.on.error=F)$value
  return(val)
}

#------------------------------------------------------------------
#FFT the interferogram and clean it up to make a pretty spectrum
#------------------------------------------------------------------
fft.interferogram<-function(interfero,Delt,dm,plot.typ=NULL) {

  leng<-2*Delt
  n<-leng/dm
  j<-seq(2,floor(n/2)+1,1)
  lambda<-leng/(j-1)
  nu.t<-1/lambda
  #FFT:
  zf<-fft(interfero)
  #Cleanup a bit:
  spectrum<-abs(zf)[2:(floor(n/2)+1)]/(n/2)

  if(plot.typ=="Absorbance") {
    #Intensity spectrum from fft of the interferogram:
    plot(nu.t/100,spectrum,typ="l",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),ylab=expression(I(tilde(nu))));
  }

  if(plot.typ=="%T") {
    #%T spectrum
    spec.idxs<-which(nu.t<=4500*100 & nu.t>=400*100)
    plot(nu.t[spec.idxs]/100,(B(nu.t[spec.idxs],1500)-spectrum[spec.idxs])/B(nu.t[spec.idxs],1500)*100,typ="h",ylab="%T",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))))
  }

  spec.idxs<-which(nu.t<=4500*100 & nu.t>=400*100)
  wns<-nu.t[spec.idxs]/100
  percTs<-(B(nu.t[spec.idxs],1500)-spectrum[spec.idxs])/B(nu.t[spec.idxs],1500)*100

  return(list(spectrum,cbind(wns,percTs)))

}

#------------------------------------------------------------------
#FFT the interferogram and clean it up to make a pretty spectrum
#------------------------------------------------------------------
fft.interferogram2<-function(interfero, plot.typ="%T") {

  #Mirror positions. Needed to compute the interferogram:
  dm<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delt<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  #x<-seq(-(Delta),Delta-dx,dx)  #delta, units: m
  #length(x) #Number of points in the interferogram


  leng<-2*Delt
  n<-leng/dm
  j<-seq(2,floor(n/2)+1,1)
  lambda<-leng/(j-1)
  nu.t<-1/lambda
  #FFT:
  zf<-fft(interfero)
  #Cleanup a bit:
  spectrum<-abs(zf)[2:(floor(n/2)+1)]/(n/2)

  if(plot.typ=="Absorbance") {
    #Intensity spectrum from fft of the interferogram:
    plot(nu.t/100,spectrum,typ="l",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),ylab=expression(I(tilde(nu))), xlim=rev(range(nu.t/100)));
  }

  if(plot.typ=="%T") {
    #%T spectrum
    spec.idxs<-which(nu.t<=4500*100 & nu.t>=400*100)
    xf <- nu.t[spec.idxs]/100
    yf <- (B(nu.t[spec.idxs],1500)-spectrum[spec.idxs])/B(nu.t[spec.idxs],1500)*100
    plot(xf, yf,typ="l",ylab="%T",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))), xlim=rev(range(xf)))
  }

  if(plot.typ=="%Th") {
    #%T spectrum
    spec.idxs<-which(nu.t<=4500*100 & nu.t>=400*100)
    xf <- nu.t[spec.idxs]/100
    yf <- (B(nu.t[spec.idxs],1500)-spectrum[spec.idxs])/B(nu.t[spec.idxs],1500)*100
    plot(xf, yf, typ="h",ylab="%T",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))), xlim=rev(range(xf)))
  }

  # I forget, why is this needed??
  spec.idxs<-which(nu.t<=4500*100 & nu.t>=400*100)
  wns<-nu.t[spec.idxs]/100
  percTs<-(B(nu.t[spec.idxs],1500)-spectrum[spec.idxs])/B(nu.t[spec.idxs],1500)*100

  spectrum <- spectrum[spec.idxs]

  #return(list(spectrum,cbind(wns,percTs)))
  spectrum.info.mat <- cbind(wns, spectrum, percTs)
  colnames(spectrum.info.mat) <- c("nu.tilde (cm-1)", "I", "%T")
  return(spectrum.info.mat)

}

#------------------------------------------
#Simulate a (random) vibrational spectrum
#------------------------------------------
simulate.a.vibrational.spectrum<-function(plotQ=TRUE, source.only=FALSE){

  #Mirror positions. Needed to compute the interferogram:
  dx<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delta<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  x<-seq(-(Delta),Delta-dx,dx)  #delta, units: m
  length(x) #Number of points in the interferogram

  #Wavenumber axis. Needed to simulate a spectrum:
  nu<-seq(from=100*400,to=100*4500-dnu,by=dnu)    #nu-tilde (abbreviated nu!), units: m^-1
  lamb<-1/nu #Convenient to invert for interferogram computation
  length(nu) #number of spectral points. NOTE: nu at Fourier freqs to avoid spectral leakage and having to window.

  #Source intensities:
  Temperature<-1500
  Bnu<-B(nu,Temperature)

  if(source.only==FALSE){
    #Simulate a spectrum by "absorbing" part of some of the intensities:
    idxs<-sample(which(nu<=4500*100 & nu>=400*100),20,replace=F) #The nu-tilde will be where the "absorbances" occur.
    Bnu[idxs]<-(Bnu[idxs]-runif(20)*Bnu[idxs])                   #This determines how much intensity is absorbed at each index

  }

  if(plotQ==TRUE){
    plot(nu/100,Bnu,typ="l",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),ylab=expression(Absorbance(tilde(nu)))) #The simulated spectrum. Hopefully after processing the output will look like this!
  }

  spectrum<-cbind(nu,Bnu)
  colnames(spectrum)<-c("nu-tilde (m^-1)","Absorbance")
  return(spectrum)
}



#------------------------------------------
#Generate a vibrational spectrum by inputting some frequencies
#------------------------------------------
generate.a.vibrational.spectrum<-function(frequency.vec, plotQ=TRUE){

  #Mirror positions. Needed to compute the interferogram:
  dx<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delta<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  x<-seq(-(Delta),Delta-dx,dx)  #delta, units: m
  length(x) #Number of points in the interferogram

  #Wavenumber axis. Needed to simulate a spectrum:
  nu<-seq(from=100*400,to=100*4500-dnu,by=dnu)    #nu-tilde (abbreviated nu!), units: m^-1
  lamb<-1/nu #Convenient to invert for interferogram computation
  length(nu) #number of spectral points. NOTE: nu at Fourier freqs to avoid spectral leakage and having to window.

  #Source intensities:
  Temperature<-1500
  Bnu<-B(nu,Temperature)


  #Simulate a spectrum by "absorbing" part of some of the intensities:
  #idxs<-NULL
  for(i in 1:length(frequency.vec)){
    #select the index of the nearest freq incrment
    idx <- which.min(abs(nu - (frequency.vec[i]*100)))
    #idxs <- c(idxs,idx)
    Bnu[idx]<-(Bnu[idx]-runif(1)*Bnu[idx])                   #This determines how much intensity is absorbed at each index
  }


  if(plotQ==T){
    plot(nu/100,Bnu,typ="l",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),ylab=expression(Absorbance(tilde(nu)))) #The simulated spectrum. Hopefully after processing the output will look like this!
  }

  spectrum<-cbind(nu,Bnu)
  colnames(spectrum)<-c("nu-tilde (m^-1)","Absorbance")
  return(spectrum)
}


#---------------------------------------------------------------------------------
#Generate an interferogram from a spectrum of freqs by adding up their sinusoids
#---------------------------------------------------------------------------------
make.interferogram<-function(spectrum, zoomQ=FALSE){

  #Mirror positions. Needed to compute the interferogram:
  dx<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delta<-1/dnu * 0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  x<-seq(-(Delta),Delta-dx,dx)  #delta, units: m
  length(x) #Number of points in the interferogram

  #Wavenumber axis. Needed to simulate a spectrum:
  #nu<-seq(from=100*400,to=100*4500-dnu,by=dnu)    #nu-tilde (abbreviated nu!), units: m^-1
  nu<-spectrum[,1]
  lamb<-1/nu #Convenient to invert for interferogram computation
  length(nu) #number of spectral points. NOTE: nu at Fourier freqs to avoid spectral leakage and having to window.

  #Construct interferogram:
  #Generate the cosine waves for each freq:
  Bnu<-spectrum[,2]
  cos.mat<-sapply(1:length(Bnu),function(i){Bnu[i]*cos(2*pi*x/lamb[i])}) #Columns are the cosine waves of intensity for each mirror position
  #Add up the cosine waves to get the interferogram:
  interferogram<-rowSums(cos.mat)
  #Raw Interferogram:
  if(zoomQ==F){
    plot(interferogram,typ="l", main="Interferogram")
  }
  if(zoomQ==T){
    #Centerburst zoomed in:
    plot(x[floor((length(x)/2)-100):floor((length(x)/2)+100)]*1000,interferogram[floor((length(x)/2)-100):floor((length(x)/2)+100)],typ="l",ylab=expression(paste("I(",delta,")")),xlab=expression(paste(delta," (mm)")),main="Interferogram")
  }

  return(interferogram)

}


#----------------------------------------------
#Mix a few sinusoids together,view interferogram
#----------------------------------------------
mix.some.freqs<-function(freq.vec, plot.typ="spectrum"){

  #Mirror positions. Needed to compute the interferogram:
  dx<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delta<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  x<-seq(-(Delta),Delta-dx,dx)  #delta, units: m
  length(x) #Number of points in the interferogram

  #Wavenumber axis. Needed to simulate a spectrum:
  nu<-seq(from=100*400,to=100*4500-dnu,by=dnu)    #nu-tilde (abbreviated nu!), units: m^-1
  lamb<-1/nu #Convenient to invert for interferogram computation
  length(nu) #number of spectral points. NOTE: nu at Fourier freqs to avoid spectral leakage and having to window.

  #Source intensities:
  Temperature<-1500
  Bnu<-B(nu,Temperature)

  #Simulate a spectrum by "absorbing" intensities at chosen wns:
  #Pick the close nus to the ones specified
  idxs<-sapply(1:length(freq.vec), function(x){max(which((nu<=(freq.vec[x]*100))))})
  #print(idxs)

  Bnu[idxs]<-(Bnu[idxs]*0)                   #This determines how much intensity is absorbed at each index
  #plot(nu/100,Bnu,typ="l",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),ylab=expression(Absorbance(tilde(nu)))) #The simulated spectrum. Hopefully after processing the output will look like this!


  spec.idxs<-which(nu<=4500*100 & nu>=400*100)
  if(plot.typ=="spectrum"){
    plot(nu[spec.idxs]/100,(B(nu[spec.idxs],1500)-Bnu[spec.idxs])/B(nu[spec.idxs],1500)*100, typ="h",ylab="%T",xlab=expression(paste(tilde(nu)," ",(cm^{-1}))))

    spectrum<-cbind(nu, Bnu)
    colnames(spectrum)<-c("nu-tilde (m^-1)", "Intensity")
    return(spectrum)
  }

  if(plot.typ=="interferogram"){
    Bnu.zeroed<-numeric(length(Bnu))
    freq.idxs<-which(Bnu==0)
    Bnu.zeroed[freq.idxs]<-1

    interferogram<-make.interferogram(cbind(nu,Bnu.zeroed),zoomQ=T)
    return(interferogram)
  }

}


#----------------------------------------------
#Plots a centerburst, zoomed in,interferogram
#----------------------------------------------
plot.interferogram<-function(interfero, zoomQ=F, zoom.clip=100){

  #CAUTION: Assumes interferogram was generated with these parameters:

  #Mirror positions. Needed to compute the interferogram:
  dx<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delta<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  xx<-seq(-(Delta),Delta-dx,dx)  #delta, units: m

  if(zoomQ==TRUE){
    plot(xx[floor((length(xx)/2)-zoom.clip):floor((length(xx)/2)+zoom.clip)]*1000,
         interfero[floor((length(xx)/2)-zoom.clip):floor((length(xx)/2)+zoom.clip)],
         typ="l",
         ylab=expression(paste("I(",delta,")")),
         xlab=expression(paste(delta," (mm)")),main="Interferogram")
  } else {
    plot(xx,
         interfero,
         typ="l",
         ylab=expression(paste("I(",delta,")")),
         xlab=expression(paste(delta," (mm)")),main="Interferogram")
  }

}


# Plot a spectrum with pretty axis labels
plot.spectrum <- function(spectrum.mat) {

  plot(spectrum.mat[,1], spectrum.mat[,2],
       typ="l",
       ylab="I",
       xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),
       xlim=rev(range(xf)))

}


# Read in and process a JDX format spectrum from the NIST WebBook
process.JDX.spectrum <- function(input.dat, printQ=FALSE, plotQ=FALSE){

  if(is.character(input.dat)){ # If a string path or URL is passed in
    in.dat <- readJDX(input.dat)
    spectrm <- cbind(in.dat[[4]]$x, in.dat[[4]]$y) # x in cm^-1
    if(printQ==TRUE){
      print(in.dat$metadata)
    }
  } else { # If a JDX spectrum is passed in two-column format
    spectrm <- input.dat
  }

  xx <- spectrm[,1] # Wavenumbers (cm^-1)
  yy <- spectrm[,2] # Absorbances

  # Pad spectrum to run from 400cm^-1 to 4500cm^-1 and Re-sample to 2cm^-1 increments:
  # Padding:
  left.padding  <- seq(from=min(xx)-ceiling((min(xx)-400)/4)*4, to=min(xx)-4, by=4)
  right.padding <- seq(from=max(xx)+4, to = max(xx) + ceiling((4500 - max(xx))/4)*4, by=4)

  xx <- c(left.padding, xx, right.padding)
  #yy <- c(numeric(length(left.padding)),yy,numeric(length(right.padding))) # Zero-pads right now. Do min y instead??
  yy <- c(rep(min(yy), length(left.padding)), yy, rep(min(yy), length(right.padding))) # Pads with the minimum Absorbance


  # Re-sampling spectrum from 4cm^-1 increments to 2cm^-1 increments:
  spec.func <- splinefun(xx, yy)
  dnu.tilde <- 2
  xx.resamp <- seq(from=400,to=4500-dnu.tilde,by=dnu.tilde)
  yy.resamp <- abs(spec.func(xx.resamp))
  spectrm.resamp  <- cbind(xx.resamp*100, yy.resamp) # put x in m^-1
  colnames(spectrm.resamp) <- c("nu-tilde (m^-1)", "Absorbance")

  if(plotQ==T){
    plot(spectrm.resamp[,1]/100,spectrm.resamp[,2],typ="l",
         xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),
         ylab=expression(Absorbance(tilde(nu)))
        )
  }

  return(spectrm.resamp)

}


# Find peaks
find.peaks<-function(spectrum.mat, spectrum.typ="Absorbance", deriv.tol, peak.tol=NULL, plotQ=FALSE) {

  # Taken from profileslib

  xx     <- spectrum.mat[,1] # Wavenumbers in cm^-1
  profil <- spectrum.mat[,2]     # Absorbances

  #fit profile with splines and compute first derivative:
  #profile "function":
  psf<-splinefun(1:length(profil),profil)

  #First derivative:
  dpsf<-psf(1:length(profil),deriv=1)

  #Find zero crossings of derivative:
  deriv.zeros<-NULL
  deriv.vals<-NULL
  for(i in 1:(length(dpsf)-1))
  {
    if(dpsf[i]*dpsf[i+1]<0)
    {
      deriv.zeros<-c(deriv.zeros,i)
    }
  }

  #Drop out the ends if they had zero derivatives:
  if(1 %in% deriv.zeros)
  {
    deriv.zeros<-deriv.zeros[-which(deriv.zeros==1)]
  }
  if(length(profil) %in% deriv.zeros)
  {
    deriv.zeros<-deriv.zeros[-which(deriv.zeros==length(profil))]
  }

  #Pluck out any small extrema. They are probably noise from the
  #differencing procedure above, or due to a long flat set of
  #adjacent extrema points.
  extrema.height.diffs<-NULL
  for(i in 2:length(deriv.zeros))
  {
    #Subtract Right (current) extremum from extremum to the Left
    val<-profil[deriv.zeros[i]]-profil[deriv.zeros[i-1]]
    extrema.height.diffs<-c(extrema.height.diffs,val)
  }

  #Get rid of "noise" extrema:
  col1<-which(abs(extrema.height.diffs)>deriv.tol)

  #Keep the first (left most) extremum:
  col1<-col1+1 #Shift the indices by 1
  col1<-c(1,col1) #Tack on index for left most extremum

  #Extrema (df/dx~0) indices:
  deriv.zeros<-deriv.zeros[col1]

  #extrema.height.differences missing difference for first extremum.
  #To the left of first extremum may be a noise extremum (e.g.profile point 1).
  #To get a more definitive idea of whether it is a min or a max, subtract
  #it from the extremum to its right
  first.extremum.height.diff<-val<-profil[deriv.zeros[1]]-profil[deriv.zeros[2]]
  extrema.height.diffs<-c(first.extremum.height.diff,extrema.height.diffs)
  extrema.height.diffs<-extrema.height.diffs[col1]

  #Max=1,Min=-1. THERE SHOULD BE NO ZEROS!!!!!!!
  extrema.typ<-sign(extrema.height.diffs)
  if(0 %in% extrema.typ)
  {
    print("Undefined max/min!")
    print(cbind(col1,extrema.height.diffs,deriv.zeros,extrema.typ))
    return(0)
  }


  if(spectrum.typ=="Absorbance"){
    max.deriv.zeros.idxs <- which(extrema.typ == 1)
  } else if(spectrum.typ=="%T"){
    max.deriv.zeros.idxs <- which(extrema.typ == -1)
  } else {
    stop("Specify spectrum.typ = Absorbance or %T!")
  }

  max.xs.idxs          <- deriv.zeros[max.deriv.zeros.idxs]
  #print(cbind(extrema.typ[max.deriv.zeros.idxs],deriv.zeros[max.deriv.zeros.idxs]))

  if(is.null(peak.tol)){
    peaks   <- xx[max.xs.idxs]
    max.yys <- profil[max.xs.idxs]
  } else {
    peaks   <- xx[max.xs.idxs]
    max.yys <- profil[max.xs.idxs]

    select.max.xs.idxs <- which(profil[max.xs.idxs] >= peak.tol)

    peaks   <- peaks[select.max.xs.idxs]
    max.yys <- max.yys[select.max.xs.idxs]
  }


  if(plotQ==TRUE) {
    ymax<-max(profil)
    ymin<-min(profil)
    plot(xx, profil, typ="l", xlim=rev(c(min(xx),max(xx))), ylim=c(ymin,ymax),
         xlab=expression(paste(tilde(nu)," ",(cm^{-1}))),
         ylab=expression(I(tilde(nu)))
        )
    #points(xx[max.xs.idxs], profil[max.xs.idxs], col="red")
    points(peaks, max.yys, col="red")
  }

  return(peaks)

}
