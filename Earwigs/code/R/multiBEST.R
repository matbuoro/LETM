
rm(list=ls())   # Clear memory

library(ggplot2)
library(mcmcplots)

##------- REPERTORY---------##
#setwd("~/Dropbox/RESEARCH/PROJECT/LETM/TEST/PEDIGREE/")
setwd('~/Dropbox/Research/PROJECT/LETM/BUZATTO/Data for Mat/EARWIGS/')

##-------------SIMULATED DATA ------------ ##

name.data=c('br','ewo','wwo', 'sta')

data.br <- new.env();load("DATA/br.data.RData",envir=data.br)
data.ewo <- new.env();load("DATA/ewo.data.RData",envir=data.ewo)
data.wwo <- new.env();load("DATA/wwo.data.RData",envir=data.wwo)
data.sta <- new.env();load("DATA/sta.data.RData",envir=data.sta)

data=list()
data[["br"]]<-data.br$data
data[["ewo"]]<-data.ewo$data
data[["wwo"]]<-data.wwo$data
data[["sta"]]<-data.sta$data

Y=list();N=NULL
for (j in 1:length(data)){
  N[j]=dim(data[[j]])[1]
  Y[[j]] = data[[j]]$Pronotum  
}

# Load the data:
y<-unlist(Y)
x<-rep(1:4,time=N)



multiBESTmcmc = function( y,x, numSavedSteps=100000, thinSteps=1){#, showMCMC=TRUE) { 
  # This function generates an MCMC sample from the posterior distribution.
  # Description of arguments:
  # showMCMC is a flag for displaying diagnostic graphs of the chains.
  #    If F (the default), no chain graphs are displayed. If T, they are.
  
  require(rjags)
  
  #------------------------------------------------------------------------------
  # THE MODEL.
  modelString = "
  model {
  for ( i in 1:Ntotal ) {
  y[i] ~ dt( mu[x[i]] , tau[x[i]] , nu )
  }
  for ( j in 1:M ) {
  mu[j] ~ dnorm( muM , muP )
  tau[j] <- 1/pow( sigma[j] , 2 )
  sigma[j] ~ dunif( sigmaLow , sigmaHigh )
  }
  nu <- nuMinusOne+1
  nuMinusOne ~ dexp(1/29)
  }
  " # close quote for modelString
  # Write out modelString to a text file
  writeLines( modelString , con="BESTmodel.txt" )
  
  #------------------------------------------------------------------------------
  # THE DATA.
  Ntotal = length(y)
  # Specify the data in a list, for later shipment to JAGS:
  dataList = list(
    y = y ,
    x = x ,
    Ntotal = Ntotal ,
    M = max(x),
    muM = mean(y) ,
    muP = 0.000001 * 1/sd(y)^2 ,
    sigmaLow = sd(y) / 1000 ,
    sigmaHigh = sd(y) * 1000 
  )
  
  #------------------------------------------------------------------------------
  # INTIALIZE THE CHAINS.
  # Initial values of MCMC chains based on data:
  mu = as.vector(aggregate( y~x, data=as.data.frame(cbind(y,x)),mean )$y)
  sigma = as.vector(aggregate( y~x, data=as.data.frame(cbind(y,x)),sd)$y)
  # Regarding initial values in next line: (1) sigma will tend to be too big if 
  # the data have outliers, and (2) nu starts at 5 as a moderate value. These
  # initial values keep the burn-in period moderate.
  initsList = list( mu = mu , sigma = sigma , nuMinusOne = 4 )
  
  #------------------------------------------------------------------------------
  # RUN THE CHAINS
  
  parameters = c( "mu" , "sigma" , "nu" )     # The parameters to be monitored
  adaptSteps = 500               # Number of steps to "tune" the samplers
  burnInSteps = 1000
  nChains = 3 
  nIter = ceiling( ( numSavedSteps * thinSteps ) / nChains )
  # Create, initialize, and adapt the model:
  jagsModel = jags.model( "BESTmodel.txt" , data=dataList , inits=initsList , 
                          n.chains=nChains , n.adapt=adaptSteps )
  # Burn-in:
  cat( "Burning in the MCMC chain...\n" )
  update( jagsModel , n.iter=burnInSteps )
  # The saved MCMC chain:
  cat( "Sampling final MCMC chain...\n" )
  codaSamples = coda.samples( jagsModel , variable.names=parameters , 
                              n.iter=nIter , thin=thinSteps )
  # resulting codaSamples object has these indices: 
  #   codaSamples[[ chainIdx ]][ stepIdx , paramIdx ]
  
  #------------------------------------------------------------------------------
  # EXAMINE THE RESULTS
#   if ( showMCMC ) {
#     openGraph(width=7,height=7)
#     autocorr.plot( codaSamples[[1]] , ask=FALSE )
#     show( gelman.diag( codaSamples ) )
#     effectiveChainLength = effectiveSize( codaSamples ) 
#     show( effectiveChainLength )
#   }
  
  # Convert coda-object codaSamples to matrix object for easier handling.
  # But note that this concatenates the different chains into one long chain.
  # Result is mcmcChain[ stepIdx , paramIdx ]
  mcmcChain = as.matrix( codaSamples )
  return( mcmcChain )
  
} # end function BESTmcmc





BESTsummary = function( M,mcmcChain,name.data) {
  source("HDIofMCMC.R")
  mcmcSummary = function( paramSampleVec , compVal=NULL ) {
    meanParam = mean( paramSampleVec )
    medianParam = median( paramSampleVec )
    dres = density( paramSampleVec )
    modeParam = dres$x[which.max(dres$y)]
    hdiLim = HDIofMCMC( paramSampleVec )
    if ( !is.null(compVal) ) {
      pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                      / length( paramSampleVec ) )
    } else {
      pcgtCompVal=NA
    }
    return( c( meanParam , medianParam , modeParam , hdiLim , pcgtCompVal ) )
  }
  # Define matrix for storing summary info:
  comb <- combn(1:M, 2)
  ncombn <- dim(comb)[2]
  
  mu=paste("mu", 1:M, sep = "")
  sigma=paste("sigma", 1:M, sep = "")
  muDiff=paste("muDiff", 1:ncombn, sep = "")
  sigmaDiff=paste("sigmaDiff", 1:ncombn, sep = "")
  effSz=paste("effSz", 1:ncombn, sep = "")
  
  summaryInfo = matrix( 0 , nrow=(length(mu) + length(sigma)+length(muDiff)+length(sigmaDiff)+length(effSz)+2) , ncol=6 , dimnames=list(
    PARAMETER=c( mu , muDiff , sigma , sigmaDiff ,"nu" , "nuLog10" , effSz ),
    SUMMARY.INFO=c( "mean" , "median" , "mode" , "HDIlow" , "HDIhigh" ,"pcgtZero" ) 
  ) )
  
  
  for (i in 1:M){
    summaryInfo[ mu[i] , ] = mcmcSummary( mcmcChain[,paste("mu[",i,"]",sep="")],compVal=0 ) 
    summaryInfo[ sigma[i] , ] = mcmcSummary( mcmcChain[,paste("sigma[",i,"]",sep="")],compVal=0 ) 
  }
  
  for (j in 1:ncombn){
  summaryInfo[ muDiff[j] , ] = mcmcSummary( mcmcChain[,paste("mu[",comb[1,j],"]",sep="")]
                                           - mcmcChain[,paste("mu[",comb[2,j],"]",sep="")] , 
                                           compVal=0 )
  
  summaryInfo[ sigmaDiff[j] , ] = mcmcSummary( mcmcChain[,paste("sigma[",comb[1,j],"]",sep="")]
                                                                - mcmcChain[,paste("sigma[",comb[2,j],"]",sep="")] , 
                                                                compVal=0 )
  
  #N1 = length(y1)
  #N2 = length(y2)
  effSzChain = ( ( mcmcChain[,paste("mu[",comb[1,j],"]",sep="")] - mcmcChain[,paste("mu[",comb[2,j],"]",sep="")] ) 
                 / sqrt( ( mcmcChain[,paste("sigma[",comb[1,j],"]",sep="")]^2 + mcmcChain[,paste("sigma[",comb[2,j],"]",sep="")]^2 ) / 2 ) ) 
  summaryInfo[ effSz[j] , ] = mcmcSummary( effSzChain , compVal=0 )
  # Or, use sample-size weighted version:
  # effSz = ( mu1 - mu2 ) / sqrt( ( sigma1^2 *(N1-1) + sigma2^2 *(N2-1) ) 
  #                               / (N1+N2-2) )
  # Be sure also to change plot label in BESTplot function, below.
  }
  

  summaryInfo[ "nu" , ] = mcmcSummary( mcmcChain[,"nu"] )
  summaryInfo[ "nuLog10" , ] = mcmcSummary( log10(mcmcChain[,"nu"]) )
 
  if(!is.null(name.data)){
    library(plyr)
    combinaison=mapvalues(comb, from = 1:M, to = name.data)    
  }


  return( list(name.data=name.data,combinaison=combinaison,summary=summaryInfo) )
}

BESTsummary(4,mcmcChain,name.data)
