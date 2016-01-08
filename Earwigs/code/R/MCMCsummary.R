
MCMCsummary = function( codaSamples) {
  
  # Result is mcmcChain[ stepIdx , paramIdx ]
  mcmcChain = as.matrix( codaSamples )
  
  # Names of parameters
  par = colnames(mcmcChain)
  
  source("HDIofMCMC.R")
  
  mcmcSummary = function( paramSampleVec , compVal=NULL ) {
    meanParam = mean( paramSampleVec )
    sdParam = sd( paramSampleVec )
    medianParam = median( paramSampleVec )
    dres = density( paramSampleVec )
    #modeParam = dres$x[which.max(dres$y)]
    hdiLim = HDIofMCMC( paramSampleVec )
    if ( !is.null(compVal) ) {
      pcgtCompVal = ( 100 * sum( paramSampleVec > compVal ) 
                      / length( paramSampleVec ) )
    } else {
      pcgtCompVal=NA
    }
    return( c( meanParam , sdParam, medianParam , hdiLim , pcgtCompVal ) )
  }
 
  # Define matrix for storing summary info:  
  summaryInfo = matrix( 0 , nrow=length(par) , ncol=6 , dimnames=list(
    PARAMETER=par,
    SUMMARY.INFO=c( "mean" ,"sd", "median" , "2.5%" , "97.5%" ,"pcgtZero" ) 
  ) )
  
  for (i in 1:length(par)){
    summaryInfo[ i , ] = mcmcSummary( mcmcChain[,par[i]],compVal=0 ) 
  }
  
  return(summaryInfo)
}