rm(list=ls())   # Clear memory

## DIRECTORY
setwd("~/Documents/RESEARCH/PROJECTS/LETM/")

## R PACKAGES
library(rjags)
library(mcmcplots)

#---------------------DATA-----------------------------------#
nf <-20    # nb batches
nind<- 20   # nb individuals /batch
n<-rep(nind,nf)
F<-rep(1:nf,each=nind)
N <- sum(n)   # Population size


h2 <- 0.25   # Heritability
Vp <- 1     # Phenotypic variance
Va <- h2 * Vp   # Additive genetic variance
Ve <- (1-h2)*Vp   # Environmental variance
SD.a <- sqrt(Va)  # Additive genetic standard deviation
SD.e <- sqrt(Ve)  # Environmental standard deviation


## 1. GENERATE DATASET

# Generate pedigree

GeneratePedigree = function(nDesc,nMother, nFather,id.Family){
  
  nAsc=nMother+nFather # Number of Ascendants
  N = nAsc + nDesc # Total number of individuals
  ped=array(NA,dim=c(N,4));colnames(ped)<-c("ID","Sire","Dam","family") # pedigree matrix
  ped[,"ID"]=1:N # index
  ped[,"family"]= id.Family # index
  id.dam=1:nMother
  id.sire=(nMother+1):(nMother+nFather)
  
  for (i in (nAsc+1):dim(ped)[1]){
    ped[i,"Sire"]=ifelse(nFather==1,id.sire,sample(id.sire,1,replace=TRUE)) 
    ped[i,"Dam"]= ifelse(nMother==1,id.dam,sample(id.dam,1,replace=TRUE))
  }
  
  return(list(nAsc=nAsc,ped=as.data.frame(ped)))
}



nfamily = length(unique(F))
pedigree=list()
for (n in 1:nfamily){
  nDesc=table(F)[n]
  #nMother=rpois(1,2)+1
  nMother=1
  #nFather=1
  nFather=rpois(1,2)+1
  nAsc=nMother + nFather
  data=GeneratePedigree(nDesc,nMother, nFather,n)
  pedigree[[n]]=data$ped
}

data <- do.call(rbind,pedigree)
data <- as.data.frame(na.omit(data))

## Re-indexing of the Sire and dam ID
data$Sire = data$Sire + data$family*100
sire.tmp <- unique(data$Sire);new.sire.id <- NULL
for (i in 1:dim(data)[1]){
  for (j in 1:length(sire.tmp)){
    if(data$Sire[i] == sire.tmp[j]) new.sire.id[i] <- j
  }
}
data$Sire <- new.sire.id

data$Dam = data$Dam + data$family*100
dam.tmp <- unique(data$Dam);new.dam.id <- NULL
for (i in 1:dim(data)[1]){
  for (j in 1:length(dam.tmp)){
    if(data$Dam[i] == dam.tmp[j]) new.dam.id[i] <- max(data$Sire) + j
  } 
}
data$Dam <- new.dam.id

nAsc=max(c(data$Dam,data$Sire)) # N parents
nDesc = dim(data)[1] # N offsprings


# THRESHOLD
a_asc <- rnorm(nAsc,0, sqrt(Va))
a=mu_a=NULL
for (i in 1:nDesc){
  mu_a[i]<-(a_asc[data$Dam[i]]+a_asc[data$Sire[i]])/2
  a[i] <-  rnorm(1, mu_a[i], sqrt(0.5*Va))
}
mu_theta <- 0  
Theta <- mu_theta + a


# ENVIRONMENTAL CUE
X_tmp <- rnorm(N,0,1)   # Observable cue
X <-as.vector(scale(X_tmp))

Eta=NULL
Eta<-rnorm(nDesc,X,sqrt(Ve)) # Proximate cue

# Calculate tactics by comparing Proximate cue to threshold
Y <- ifelse(Eta > Theta, 1, 0)


data.jags<-list(
  X = X,
  Y = Y,
  Ndesc = length(Y),    # Population size
  Nasc = nAsc,
  Dam = data$Dam, # ID of dams
  Sire = data$Sire # ID of sires
)

#---------------------MODEL-----------------------------------#

sink("LETM.md")
cat("

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
  #tau_T~dchisqr(1)
  #sigma2_T <- 1/(tau_T/(sd(X)*sd(X)))
  
  # variance of the latent threshold
  sigmaT_theta <- sigma2_T * h2
  
  # variance of the proximate cue
  sigma2_eta <- sigma2_T * (1-h2)
  
  # Heritability
  h2~dbeta(1,1)
  
}        # End model
    ",fill = TRUE)
sink()


#---------------------ANALYSIS-----------------------------------#

nchains = 2 # Number of chains to run.
nadapt= 1000
nburnin=5000
nstore=10000 # Total number of steps in chains to save.
nthin=1 # Number of steps to "thin" (1=keep every step).

## PARAMETERS
parameters<-c("mu_theta","sigma2_T","sigma2_eta","sigma2_theta","h2")


## INITIAL VALUES
mu_theta.inits=mean(X)
eta.inits=ifelse(Y==0,runif(1,mu_theta-0.1,mu_theta-0.01),runif(1,mu_theta+0.01,mu_theta+0.1))

inits = function(){list(
  mu_theta=mu_theta.inits,
  eta=eta.inits,
  #tau_T=runif(1,0,0.1),
  sigma2_T=runif(1,0,10),
  h2=rbeta(1,2,2)
)}

burnin = function(iter){
  if (iter<=10000) {nburnin = iter * 0.2} else { nburnin = 2500};
  return(nburnin)
}


## ANALYSIS (Using Jags)
cat( "Adapt...\n" )
model.fit <- jags.model( file = "LETM.md", 
                         data = data.jags,
                         inits,
                         n.chains = nchains,n.adapt = nadapt )

cat( "Burnin...\n" )
update(model.fit, n.iter=burnin(nstore))

cat( "Sampling MCMC....\n" )
model.mcmc <- coda.samples(model.fit,
                           var = parameters,
                           n.iter = nstore*nthin,thin = nthin)

#---------------------RESULTS-----------------------------------#
summary(model.mcmc)

traplot(model.mcmc, "h2");

caterplot(model.mcmc, "h2"); caterpoints(h2)
caterplot(model.mcmc, "mu_theta"); caterpoints(mu_theta)
caterplot(model.mcmc, "sigma2_theta"); caterpoints(Va)
caterplot(model.mcmc, "sigma2_eta"); caterpoints(Ve)

