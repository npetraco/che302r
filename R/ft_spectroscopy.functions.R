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
#FFT the interferogram and clean it up to make a prettry spectrum
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
  Delta<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
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
plot.interferogram<-function(interfero){

  #CAUTION: Assumes interferogram was generated with these parameters:

  #Mirror positions. Needed to compute the interferogram:
  dx<-0.000001 #d-delta units: m
  dnu<-2*100   #d-nu-tilde, units: m^-1
  Delta<-1/dnu *0.5 #Max discplacement of mirror. Make less than theoretical limit to keep-out aliasing down stream in the interferogram and resulting spectrum.
  x<-seq(-(Delta),Delta-dx,dx)  #delta, units: m

  plot(x[floor((length(x)/2)-100):floor((length(x)/2)+100)]*1000,
       interfero[floor((length(x)/2)-100):floor((length(x)/2)+100)],
       typ="l",ylab=expression(paste("I(",delta,")")),
       xlab=expression(paste(delta," (mm)")),main="Interferogram")
}

