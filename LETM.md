

model {
  #======== LIKELIHOOD ========#
  
  for (i in 1:Ndesc) {
    
    # Phenotypes
    Y[i]~dinterval(eta[i], theta[i]) # works with JAGS only!
    #Y[i]~dbern(p[i])
    #p[i]<-phi(z[i])
    #z[i]<-((X[i] - theta[i]) / sqrt(sigma2_eta))
    
    ## Threshold
    theta[i]<-mu_theta  + a[i]
    
    # Additive genetic value:
    a[i]~dnorm( mu_a[i] , tau_a)
    
    # Mean of the additive genetic value:
    mu_a[i] <- (a_asc[Dam[i]] + a_asc[Sire[i]]) / 2 
    
    ## Proximate cue (eta) / Environmental cue (X):
    eta[i]~dnorm(X[i], 1/sigma2_eta)
    
  } # End of loop i
  
  # Precision of the additive genetic value:
  tau_a <- (1/ (0.5*sigma2_theta)) # for offsprings
  for (j in 1:Nasc) { a_asc[j]~dnorm(0,tau_theta) } # parents
  tau_theta <- 1/sigma2_theta
  
  # Mean and sd of the proximate cue
  mu_eta <- mean(eta)
  sd_eta <- sd(eta)
  
  #======== PRIORS ========#    
  # Mean of the threshold
  mu_theta~dnorm(0, 0.001)
  
  # Total variance (sigma2_T= sigma2_theta + sigma2_eta)
  sigma2_T~dunif(0,100)
  # alternative prior:
  #tau_T~dchisqr(1)
  #sigma2_T <- 1/(tau_T/(sd(X)*sd(X)))
  
  # variance of the latent threshold
  sigmaT_theta <- sigma2_T * h2
  
  # variance of the proximate cue
  sigma2_eta <- sigma2_T * (1-h2)
  
  # Heritability
  h2~dbeta(1,1)
  
}        # End model
    
