####################################################################
# Investigating the genetic architecture of conditional strategies using the Environmental Threshold Model
#
# Bruno A. Buzatto*,1,5, Mathieu Buoro*2,3, Wade N. Hazel4 and Joseph L. Tomkins1,6
# 1 Centre for Evolutionary Biology, School of Animal Biology (M092), The University of Western Australia, 35 Stirling Highway, Crawley 6009, WA, Australia
# 2 INRA, UMR 1224, Ecologie Comportementale et Biologie des Populations de Poissons, Saint Pêe sur Nivelle, France
# 3 Univ Pau & Pays Adour, UMR 1224, Ecologie Comportementale et Biologie des Populations de Poissons, UFR Sciences et Techniques de la Côte Basque, Anglet, France
# 4 Department of Biology, DePauw University, Greencastle, IN 46135, USA 
# 5 E-mail: bruno.buzatto@gmail.com
# * These authors contributed equally to this work.
####################################################################

rm(list=ls())   # Clear memory



##------- REPERTORY---------##
path<-c('~/Documents/RESEARCH/PROJECTS/LETM/Earwigs/')
setwd(path)

##------- PACKAGES---------##
library(rjags)
library(mcmcplots)

##------------- DATA ------------ ##
#load("ISLANDS.RData") # contains data for 3 populations
data<-list()
data[["wwo"]] <- read.csv("data/DATA_ISLANDS_wwo.csv",sep=";")
data[["br"]] <- read.csv("data/DATA_ISLANDS_br.csv",sep=";")
data[["ewo"]] <- read.csv("data/DATA_ISLANDS_ewo.csv",sep=";")
population=c('wwo','br','ewo')


N=NULL;nsire=ndam=NULL
for (j in 1:length(population)){
  N[j]=dim(data[[j]])[1] # nb of individuals / population
  nsire[j]=length(unique(data[[j]]$Sire)) # nb of males / population
  ndam[j]=length(unique(data[[j]]$Dam)) # nb of females / population
}

## combining all datasets in one
data_all.tmp = do.call(rbind,data)
data_all <- cbind(data_all.tmp,I=rep(1:length(population),time=N))


## Re-index of the Sire and dam ID for the whole dataset
# here, we add j (index of popualtion: 1 to 3) *10000 to make sure all sire and dams of each populations have different ID
for (j in 1:length(population)){
  data_all$ID <- ifelse(data_all$I==j, data_all$ID + (j*10000),data_all$ID)
  data_all$Sire <- ifelse(data_all$I==j, data_all$Sire + (j*10000),data_all$Sire)
  data_all$Dam <- ifelse(data_all$I==j, data_all$Dam + (j*10000),data_all$Dam)
}

sire.tmp <- unique(data_all$Sire);new.sire.id <- NULL
dam.tmp <- unique(data_all$Dam);new.dam.id <- NULL

for (i in 1:dim(data_all)[1]){
  for (j in 1:length(sire.tmp)){
    if(data_all$Sire[i] == sire.tmp[j]) new.sire.id[i] <- j
  }
  
  for (j in 1:length(dam.tmp)){
    if(data_all$Dam[i] == dam.tmp[j]) new.dam.id[i] <- length(sire.tmp) + j
  } 
}
data_all$Sire <- new.sire.id
data_all$Dam <- new.dam.id


#### 
Y=data_all$Morph # alternative phenotypes 1/0
X=data_all$Pronotum # Observable cue
I=data_all$I # Indicator of populations

data.jags.all <- list(
  npop=length(population),
  I=I,
  n=table(data_all$I),
  J=rep(1:length(population),
  time=(nsire+ndam)),
  X=X,
  Y=Y,
  Ndesc=dim(data_all)[1], 
  Nasc=max(data_all$Sire,data_all$Dam),
  Sire=as.vector(data_all$Sire),
  Dam=as.vector(data_all$Dam)
  )


#---------------------------LETM MODEL----------------------------------#
write("
      
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
      h2[j]~dbeta(2,2)
      
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
      ", "code/models/LETM_ISLANDS.md")


#---------------------------ANALYSIS----------------------------------#

## PARAMETERS TO SAVE
parameters<-c(
              "mu_theta","sigma2_T","sigma2_eta","sigma2_theta","h2","mu_eta","sd_eta",
              "Diff_mu_theta","Diff_sigma2_theta","Diff_mu_eta","Diff_sigma2_eta","Diff_sigma2_T","Diff_h2"
              )

nchains = 2 # Number of chains to run.
nadapt= 1000
nburnin=5000
nstore=10000 # Total number of steps in chains to save.
nthin=1 # Number of steps to "thin" (1=keep every step).

#===========================
burnin = function(iter){
  if (iter<=10000) {nburnin = iter * 0.2} else { nburnin = 2500};
  return(nburnin)
}

## INITS
mu_theta=mu_theta=c(mean(X[I==1]),mean(X[I==2]),mean(X[I==3]))
eta=NULL
# simualting individual proximate cues
for (i in 1:dim(data_all)[1]){
eta[i]=ifelse(data.jags.all$Y[i]==0,runif(1,mu_theta[data.jags.all$I[i]]-0.1,mu_theta[data.jags.all$I[i]]-0.01),runif(1,mu_theta[data.jags.all$I[i]]+0.01,mu_theta[data.jags.all$I[i]]+0.1))
}

inits<-function(){list(
  eta=eta,
  mu_theta=mu_theta,
  #tau_T=runif(3,0,0.1),
  sigma2_T=runif(3,0,10),
  h2=rbeta(3,2,2)
)}



cat( "Adapt...\n" )
model.fit <- jags.model( file = "code/models/LETM_ISLANDS.md", 
                         data = data.jags.all,
                         inits,
                         n.chains = nchains,n.adapt = nadapt )

cat( "Burnin...\n" )
update(model.fit, n.iter=burnin(nstore))

# Save current states
#inits <- model.fit$state(internal=TRUE)

cat( "Sampling MCMC....\n" )
codaSamples <- coda.samples(model.fit,
                                var = parameters,
                                n.iter = nstore*nthin,thin = nthin)

save(codaSamples,file=paste('results/RESULTS_LETM_ISLANDS.RData',sep=""))


# To save individuals thresholds and cues:
# codaSamples.states <- coda.samples(model.fit,
#                             var = c("theta","eta"),
#                             n.iter = 1000*10,thin = 10)
# 
# save(codaSamples.states,file=paste('results/RESULTS_LETM_ISLANDS_States.RData',sep=""))


#==========Outputs=================
# summary(codaSamples)
# library(mcmcplots)
# traplot(codaSamples, "h2")
# caterplot(codaSamples, "h2")



