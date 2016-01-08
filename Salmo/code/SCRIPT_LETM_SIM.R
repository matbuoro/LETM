#################################################################################
### Assessing adaptive phenotypic plasticity by means of conditional          ###
### strategies from empirical data: the latent environmental threshold model  ###
###     							      ###
### Mathieu Buoro 1,2,4,  Olivier Gimenez 1,6 and Etienne Prévost 2,3,5       ###
### 									      ###
### 1 Centre d’Ecologie Fonctionnelle et Evolutive, campus CNRS,	      ###
### UMR 5175, 1919 Route de Mende, 34293 Montpellier Cedex 5, France.         ###
### 2 INRA, UMR Ecobiop, Quartier Ibarron 64310 Saint Pée s/ Nivelle, France  ###
### 3 Université de Pau et Pays de l'Adour, UMR Ecobiop, Campus de Montaury,  ###
###  64600 Anglet, France						      ###
### 4 E-mail: mathieu.buoro@cefe.cnrs.fr				      ###
### 5 E-mail: eprevost@st-pee.inra.fr					      ###
### 6 E-mail: olivier.gimenez@cefe.cnrs.fr				      ###
#################################################################################

##-----------UPDATE!!!!!--------------##
## NOW USING JAGS http://mcmc-jags.sourceforge.net/ ###

rm(list=ls())   # Clear memory

##------- REPERTORY---------##
setwd("~/Documents/RESEARCH/PROJECTS/LETM/Salmo/")


##------- LIBRARY---------##
require(MASS)
require(rjags)
require(R2jags) 
#source("http://bioconductor.org/biocLite.R")
#biocLite("GeneticsPed")
#require(GeneticsPed)
#require(magic)

#---------------------------LETM MODEL----------------------------------#
write("
      model {
      #======== LIKELIHOOD ========#
      
      for (i in 1:N) {
      
      # Phenotypes
      Y[i]~dbern(p[i])
      p[i]<-phi(z[i])
      z[i]<-((X[i] - theta[i]) / sqrt(sigma2_eta))
      
      # Threshold
      theta[i]<-mu_theta + a[i]
      
      } # End of loop i
      
      # Additive genetic value:
      a[1:N]~dmnorm(mu_a[], invA[,] * (1/sigma2_theta))
      
      
      # Using Cholesky Decomposition:
      #       a[1:N]~dmnorm(mu_a[], invA_chol[,] * (1/sigma2_theta))
      #       for ( i in 1: N){
      #       for ( j in 1: N){
      #       A_chol[i,j] <- inprod(L.theta[i,1:N], L.theta[j,1:N])
      #       }}
      #       invA_chol[1:N, 1:N] <- inverse(A_chol[,])
      
      #======== PRIORS ========#
      
      # Mean of the threshold
      mu_theta~dnorm(0, 0.001)
      
      # Mean of the additive genetic value:
      for (i in 1:N) { mu_a[i]<-0 }
      
      # Total variance (sigma2_T= sigma2_theta + sigma2_eta)
      sigma2_P~dunif(0,1000)
      
      # variance of the latent threshold
      sigma2_theta <- sigma2_P * h2
      
      # variance of the proximate cue
      sigma2_eta <- sigma2_P * (1-h2)
      
      # Heritability
      h2~dbeta(1,1)
      
      }        # End model
      ", "code/models/LETM.txt")




##-------------SIMULATED DATA ------------ ##

nSim=1 # Number of simulated data

#set.seed(10)		# Initiate random number generator
nf <-20			# nb batches
nind<- 20		# nb individuals /batch

n<-rep(nind,nf)
K<-rep(1:nf,each=nind)
N <- sum(n)		# Population size

h2 <- 0.5		# Heritability
Vp <- 1			# Phenotypic variance
Va <- h2 * Vp		# Additive genetic variance
Ve <- (1-h2)*Vp		# Environmental variance
SD.a <- sqrt(Va)	# Additive genetic standard deviation
SD.e <- sqrt(Ve)	# Environmental standard deviation
mu <- 0			# Mean genetic value


### 1. Creating genetic relationship matrix A ###

#################################
## Fullsibs / Halfsibs Designs ##
#################################
## Model choice
model.name<-'Fullsibs'
#model.name<-'Halfsibs'

fullsibs<-0.5; halfsibs <- 0.25 	# Relatedness (0.5 for parent offspring pairs and full-siblings, 0.25 for half siblings, 0.125 for first cousins etc).
Relatedness <-ifelse(model.name=='Fullsibs',fullsibs,halfsibs)

A=array(0,dim=c(N,N))
for (i in 1:N){
  for (j in 1:N){
    ifelse(K[i]==K[j],A[i,j]<-fullsibs,0)
  }
}
diag(A)<-1

# Additive genetic variance matrix
CovA <- matrix(0,N,N)
CovA <- A * Va	

Theta <- array(,dim=c(N,nSim))
X <- array(,dim=c(N,nSim))
Eta <- array(,dim=c(N,nSim))
Y <- array(,dim=c(N,nSim))

data<-list()
for (d in 1:nSim){ # loop for each simulated datasets			
  
  # THRESHOLD (Additive genetic values / breeding values)
  Theta[,d] <- mvrnorm(1, mu=rep(0,N), Sigma=CovA)
  
  # ENVIRONMENTAL CUE
  X[,d] <- rnorm(N,0,1); X[,d]<-scale(X[,d]) # Observable cue
  Eta[,d]<-rnorm(N,X[,d],SD.e)		# Proximate cue
  
  # Determined tactics by comparing Proximate cue to threshold
  Y[,d] <- ifelse(Eta[,d] > Theta[,d], 1, 0)
  
  
  ## SAVE DATA
  data[[d]]<-list("N"=N,"Y"=Y[,d],"X"=round(X[,d],2),
             #L.theta=as.matrix(chol(A)) ## Using decomposition of Cholesky
             'invA'=as.matrix(solve(A)) # inverse of matrix A
  )
  
} #end loop d


## -> Go to "Analysis" section



####################
## Mixture Design ##
####################

model.name= "Mixture"

## As the pedigree of the fish sampled is unknown, we generated 20 mixtures’ genetic structure (i.e., additive genetic relationship matrix A) 
#source("http://bioconductor.org/biocLite.R")
#biocLite("GeneticsPed")
#require(GeneticsPed)
#require(magic)

nId=10 # Number of descendant for each batches

nFather=matrix(,nf,nSim)
nMother=matrix(,nf,nSim)
Theta <- array(,dim=c(N,nSim))
X <- array(,dim=c(N,nSim))
Eta <- array(,dim=c(N,nSim))
Y <- array(,dim=c(N,nSim))


data<-list()
for (d in 1:nSim){ # loop for each simulated datasets
  
  ## MIXTURE
  ped<-list(NULL)
  I<-list(NULL);mat<-list(NULL)
  ##invA<-list()
  for (k in 1:nf){ # loop for each batches k (1 to length(nind))
    ## sampling of number of father and mother for each batches:
    ## Note that the number of potential fathers must be higher than number of females
    repeat{nFather[k,d]=rpois(1,2);if(nFather[k,d]>0){break}}
    repeat{nMother[k,d]=rpois(1,2);if(nMother[k,d]>0 & nMother[k,d]<=nFather[k,d] & nMother[k,d]+nFather[k,d]<=n[k]){break}}
    ## Generating pedigree for each batches k using "generatePedigree" function in R packages "Pedigree"  
    ped[[k]]<-generatePedigree(nId=n[k],nGeneration=2,nFather=nFather[k,d],nMother=nMother[k,d])
    I[[k]]<-relationshipAdditive(ped[[k]]) ## Creating relationship matrix from the pedigree (ped)
    mat[[k]]<-I[[k]][-c(1:n[k]),-c(1:n[k])] ## We excluded ascendant from the matrix I
    ##invA[[k]]<-solve(A[[k]])
  } # end loop k
  
  
  # Re-Creating genetic relationship matrix A with each matrix mat in diagonal and 0 otherwise
  A<-NULL
  A<-adiag(mat[[1]],mat[[2]],mat[[3]],mat[[4]],mat[[5]],mat[[6]],mat[[7]],mat[[8]],mat[[9]],mat[[10]],
           mat[[11]],mat[[12]],mat[[13]],mat[[14]],mat[[15]],mat[[16]],mat[[17]],mat[[18]],mat[[19]],mat[[20]])
  
  
  # Additive genetic variance matrix
  CovA <- matrix(0,N,N)
  CovA <- A * Va				
  
  # THRESHOLD (Additive genetic values / breeding values)
  Theta[,d] <- mvrnorm(1, mu=rep(0,N), Sigma=CovA)
  
  # ENVIRONMENTAL CUE
  X[,d] <- rnorm(N,0,1); #X[,d]<-scale(X_tmp[,d]) # Observable cue
  Eta[,d]<-rnorm(N,X[,d],SD.e)		# Proximate cue
  
  # Determined tactics by comparing Proximate cue to threshold
  Y[,d] <- ifelse(Eta[,d] > Theta[,d], 1, 0)
  
  
  ## SAVE DATA
  data[[d]]<-list("N"=N,"Y"=Y[,d],"X"=round(X[,d],2),
             #L.theta=as.matrix(chol(A)) ## Using decomposition of Cholesky
             'invA'=as.matrix(solve(A)) # inverse of matrix A
  )
  
} #end loop d

## -> Go to "Analysis" section







#---------------------------ANALYSIS----------------------------------#

## PARAMETERS TO SAVE
parameters<-c("mu_theta","sigma2_P","sigma2_eta","sigma2_theta","h2")

nChains = 2 # Number of chains to run.
adaptSteps = 1000 # Number of steps to "tune" the samplers.
burnInSteps = 5000 # Number of steps to "burn-in" the samplers.
numSavedSteps=20000 # Total number of steps in chains to save.
thinSteps=1 # Number of steps to "thin" (1=keep every step).
nPerChain = ceiling( ( numSavedSteps * thinSteps ) / nChains ) # Steps per chain.


## INITS
inits<-function(){list(
  mu_theta=mean(X),
  sigma2_P=runif(1,0,1000),
  h2=rbeta(1,2,2)
)}


#################
#### ANALYSIS ###
#################
model.fit<-list(); fit.mcmc<-list()
for (d in 1:nSim){ # loop for each simulated datasets in MIXTURE design only!!!
  
  # ### JAGS ####
  ### Start of the run ###
  cat( "Sampling MCMC chain using R2JAGS...\n" )
  model.fit[[d]]<- jags(
    model.file = "code/models/LETM.txt",
    data[[d]],inits,parameters,
    n.chains = nChains, n.iter = nPerChain*nChains, n.burnin = burnInSteps,n.thin=thinSteps,
    progress.bar = "text")
  fit.mcmc[[d]] <- as.mcmc(model.fit[[d]])
  
} # end loop d

save(data,fit.mcmc,file=paste('results/RESULTS_LETM_SIM_',model.name,'.RData',sep=""))
