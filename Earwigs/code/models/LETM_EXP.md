
  
      model {
      #======== LIKELIHOOD ========#
      
      for (i in 1:Ndesc) {
      
      # Phenotypes
      Y[i]~dinterval(eta[i], theta[i])
      
      ## Threshold
      theta[i]<-mu_theta + a[i]
      
      # Additive genetic value of descendants:
      a[i]~dnorm( mu_a[i] , tau_a)
      
      # Mean of the additive genetic value:
      mu_a[i] <- (a_asc[Dam[i]] + a_asc[Sire[i]]) / 2 
      
      ## Proximate cue (eta) / Environmental cue (X):
      eta[i]~dnorm(X[i], 1/sigma2_eta[diet[i]])
      
      } # End of loop i
      
      tau_a <- (1/ (0.5*sigma2_THETA))
      
      # Additive genetic value of ascendants:
      for (i in 1:Nasc) { a_asc[i]~dnorm(0,1/sigma2_THETA) }
      
  
      #======== PRIORS ========#  

# Same threshold distribution whatever the diet treatments
      # Mean of the threshold
      mu_theta~dnorm(0, 0.001) 
      
      # variance of the latent threshold
      sigma2_THETA <- sigma2_theta[1]   

# Different proximate cue distributions for each the diet treatments
      for (j in 1:ndiet) {
      # variance of the proximate cue
      sigma2_eta[j] <- sigma2_T[j] * (1-h2[j])
      
      # Heritability
      h2[j]~dbeta(2,2)
      }   

      # Total variance (sigma2_T= sigma2_theta + sigma2_eta)
      #tau_T~dchisqr(1) # alternative prior
      #sigma2_T[1] <- 1/(tau_T/(sd(X)*sd(X)))
      
      sigma2_T[1]~dunif(0,10) 
      sigma2_theta[1] <- sigma2_T[1] * h2[1]

      sigma2_T[2] <- sigma2_theta[2] * (1 / h2[2])
      sigma2_theta[2]<- sigma2_theta[1] 
      
### Comparing distribution parameters between diet treatments   
      mu_eta[1] <- mean(eta[1:n[1]]) # Treatment high
      mu_eta[2] <- mean(eta[(n[1]+1):Ndesc]) # Treatment low
      mu_eta[3] <- mean(eta[])
      
      sd_eta[1] <- sd(eta[1:n[1]])  # Treatment high
      sd_eta[2] <- sd(eta[(n[1]+1):Ndesc])  # Treatment low
      sd_eta[3] <- sd(eta[])

      Diff_mu_eta <- mu_eta[1] - mu_eta[2] # Treatment high -low
      Diff_sigma2_eta <- sigma2_eta[1] - sigma2_eta[2]
      # Diff_sigma2_T <- sigma2_T[1] - sigma2_T[2]
      # Diff_h2 <- h2[1] - h2[2]
      
      }        # End model
      
      
