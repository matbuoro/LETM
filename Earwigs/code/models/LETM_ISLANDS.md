
      
      model {
      #======== LIKELIHOOD ========#
      
      for (i in 1:Ndesc) {
      
      # Phenotypes
      Y[i]~dinterval(eta[i], theta[i])
      
      ## Threshold
      theta[i]<-mu_theta[I[i]]  + a[i]

      # Additive genetic value:
      a[i]~dnorm( mu_a[i] , tau_a[I[i]])
      
      # Mean of the additive genetic value:
      mu_a[i] <- (a_asc[Dam[i]] + a_asc[Sire[i]]) / 2 
      
      ## Proximate cue (eta) / Environmental cue (X):
      eta[i]~dnorm(X[i], 1/sigma2_eta[I[i]])
  
      } # End of loop i

      # Precision of additive genetic value:      
      for (j in 1:npop){
      tau_a[j] <- (1/ (0.5*sigma2_theta[j]))
      tau_theta[j] <- 1/sigma2_theta[j]
      }
      
      # Additive genetic value of ascendants: 
      for (i in 1:Nasc) { a_asc[i]~dnorm(0,tau_theta[J[i]]) }
      
      # Calculate mean and sd for proximate cue
      mu_eta[1] <- mean(eta[1:n[1]]) # wwo
      mu_eta[2] <- mean(eta[(n[1]+1) : sum(n[1:2])]) # br
      mu_eta[3] <- mean(eta[(sum(n[1:2])+1) : sum(n[1:3])]) # ewo

      sd_eta[1] <- sd(eta[1:n[1]])
      sd_eta[2] <- sd(eta[(n[1]+1) : sum(n[1:2])])
      sd_eta[3] <- sd(eta[(sum(n[1:2])+1) : sum(n[1:3])])

      #======== PRIORS ========#  

      # Total Phenotypic variance
      #sigma2_T[1] <- 1/(tau_T[1]/(sd(X[1:n[1]])*sd(X[1:n[1]])))
      #sigma2_T[2] <- 1/(tau_T[2]/(sd(X[(n[1]+1) : sum(n[1:2])])*sd(X[(n[1]+1) : sum(n[1:2])])))
      #sigma2_T[3] <- 1/(tau_T[3]/(sd(X[(sum(n[1:2])+1) : sum(n[1:3])])*sd(X[(sum(n[1:2])+1) : sum(n[1:3])])))

for (j in 1:npop){

      # Mean of the threshold
      mu_theta[j]~dnorm(0, 0.001)
      
      # Total variance (sigma2_T= sigma2_theta + sigma2_eta)
      sigma2_T[j]~dunif(0,100) # alternative prior
      #tau_T[j]~dchisqr(1)
      
      # variance of the latent threshold
      sigma2_theta[j] <- sigma2_T[j] * h2[j]
      
      # variance of the proximate cue
      sigma2_eta[j] <- sigma2_T[j] * (1-h2[j])
      
      # Heritability
      h2[j]~dbeta(1.01,1.01)
      
}

### Comparing distribution parameters between populations
Diff_mu_theta[1] <- mu_theta[1] - mu_theta[2] # wwo vs br
Diff_mu_theta[2] <- mu_theta[1] - mu_theta[3] # wwo vs ewo
Diff_mu_theta[3] <- mu_theta[2] - mu_theta[3] # br vs ewo

Diff_sigma2_theta[1] <- sigma2_theta[1] - sigma2_theta[2] # wwo vs br
Diff_sigma2_theta[2] <- sigma2_theta[1] - sigma2_theta[3] # wwo vs ewo
Diff_sigma2_theta[3] <- sigma2_theta[2] - sigma2_theta[3] # br vs ewo

Diff_mu_eta[1] <- mu_eta[1] - mu_eta[2] # wwo vs br
Diff_mu_eta[2] <- mu_eta[1] - mu_eta[3] # wwo vs ewo
Diff_mu_eta[3] <- mu_eta[2] - mu_eta[3] # br vs ewo

Diff_sigma2_eta[1] <- sigma2_eta[1] - sigma2_eta[2]
Diff_sigma2_eta[2] <- sigma2_eta[1] - sigma2_eta[3]
Diff_sigma2_eta[3] <- sigma2_eta[2] - sigma2_eta[3]

Diff_h2[1] <- h2[1] - h2[2]
Diff_h2[2] <- h2[1] - h2[3]
Diff_h2[3] <- h2[2] - h2[3]

Diff_sigma2_T[1] <- sigma2_T[1] - sigma2_T[2]
Diff_sigma2_T[2] <- sigma2_T[1] - sigma2_T[3]
Diff_sigma2_T[3] <- sigma2_T[2] - sigma2_T[3]


      }        # End model
      
