#################################################################################
### Assessing adaptive phenotypic plasticity by means of conditional          ###
### strategies from empirical data: the latent environmental threshold model  ###
### 									      ###
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
      ", "LETM.txt")




#################################
## Fullsibs / Fullsibs Designs ##
#################################
## Model choice
model.name<-'Fullsibs'
#model.name<-'Halfsibs'

##----------------- DATA-----------------##
dat <- read.csv("data/DATA_SAT_06.csv",h=T,sep=";")
N <- length(dat$Migrant) # Population size
Y <- dat$Migrant # Alternative phenotype
X <- dat$Size # Observable cue 
ped <- dat$Site # Family

### 1. Creating genetic relationship matrix A ###
fullsibs<-0.5; halfsibs <- 0.25 	# Relatedness (0.5 for parent offspring pairs and full-siblings, 0.25 for half siblings, 0.125 for first cousins etc).
Relatedness <-ifelse(model.name=='Fullsibs',fullsibs,halfsibs)

A=array(0,dim=c(N,N))
 for (i in 1:N){
 for (j in 1:N){
 ifelse(ped[i]==ped[j],A[i,j]<-Relatedness,0)
 }
 }
diag(A)<-1


## SAVE DATA
data<-list('N'=N,'Y'=Y,'X'=X,
           #L.theta=as.matrix(chol(A)) ## Using decomposition of Cholesky
           'invA'=as.matrix(solve(A)) # inverse of matrix A
)




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



### Start of the run ###
start.time = Sys.time(); cat("Start of the run\n"); 

# ### JAGS ####
require(R2jags) 
### Start of the run ###
cat( "Sampling MCMC chain using R2JAGS...\n" )
model.fit<- jags(
  model.file = "code/models/LETM.txt",
  data,inits,parameters,
  n.chains = nChains, n.iter = nPerChain*nChains, n.burnin = burnInSteps,n.thin=thinSteps,
  progress.bar = "text")
fit.mcmc <- as.mcmc(model.fit)

# duration of the run 
end.time = Sys.time()
elapsed.time = difftime(end.time, start.time, units='mins')
cat("Sample analyzed after ", elapsed.time, ' minutes\n')

save(elapsed.time,data,model.fit,fit.mcmc,
     file=paste("results/RESULTS_LETM_SAT_",model.name, # Name of the model
       '.RData',sep=""))

#---------------------------RESULTS----------------------------------#

# # or to use some plots in coda
# # use as.mcmmc to convert rjags object into mcmc.list
# require(coda)
# summary(fit.mcmc)
# 
# # # Convert from JAGS output (coda) to ggmcmc format
# require(ggmcmc)
# # Now the ``fit.mcmc'' object (which is a coda object) must be converted into a data frame with certain attributes (number of chains, parameters, iterations, and whether parallel computing can be used) that the graphical ggmcmc functions will understand. The ggs() function does that.
# D <- ggs(fit.mcmc, parallel = FALSE)
# ggmcmc(D, file = paste("RESULTS_LETM_SAT_",model.name,".pdf",sep=""), param.page = 4)








####################
## Mixture Design ##
####################

model.name= "Mixture"

## As the pedigree of the fish sampled is unknown, we generated 20 mixtures’ genetic structure (i.e., additive genetic relationship matrix A) 
source("http://bioconductor.org/biocLite.R")
biocLite("GeneticsPed")
require(GeneticsPed)
require(magic)

#nf<-length(unique(ped)) # Number of batches
nf<-length(unique(na.omit(as.vector(ped)))) # Nb of batches
nind=NULL
for (k in 1:max(ped)){nind[k]<-length(Y[ped==k])} # Number of individuals / batches
nind<-nind[nind!=0]

nSim <-20 # Number of simulated data
nId=10
mat<-list(NULL); data<-list()

nFather=matrix(,length(nind),nSim) # Number of potential father
nMother=matrix(,length(nind),nSim)

pedigree<-array(,dim=c(nId*2,5,length(nind)))
I<-array(,dim=c(nId*2,nId*2,length(nind)))
J<-array(,dim=c(nId,nId,length(nind)))


for (d in 1:nSim){ # loop for each simulated datasets

for (k in 1:length(nind)){ # loop for each batches k (1 to length(nind))

## sampling of number of father and mother for each batches:
## Note that the number of potential fathers must be higher than number of females
repeat{nFather[k,d]=rpois(1,2);if(nFather[k,d]>0){break}}
repeat{nMother[k,d]=rpois(1,2);if(nMother[k,d]>0 & nMother[k,d]<=nFather[k,d]){break}}

## Generating pedigree for each batches k using "generatePedigree" function in R packages "Pedigree"
nId=10 # Number of descendant for each batches
pedigree[,,k]<-as.matrix(generatePedigree(nId=nId,nGeneration=2,nFather=nFather[k,d],nMother=nMother[k,d]))
colnames(pedigree)<-c('id','father','mother','generation','sex')
rownames(pedigree)<-seq(1:c(2*nId))
ped_tmp=NULL
ped_tmp<-as.Pedigree(as.matrix(pedigree[,,k]), subject="id", ascendant=c("father", "mother"),generation='generation')
ped_tmp<-sort(ped_tmp, by="generation")

## Creating relationship matrix from the pedigree (ped_tmp)
I[,,k]<-relationshipAdditive(ped_tmp)
## We excluded ascendant from the matrix I
J[,,k]<-I[-c(1:nId),-c(1:nId),k]

# We randomly pick up the number of descendant in each batches in accordance with the number of juveniles observed (nind)
smpl=NULL
smpl<-c(sample(c(1:nId),nind[k]))
mat[[k]]<-as.matrix(J[smpl,smpl,k]) # New relationship matrix for each batches

} # end of loop k

# Re-Creating genetic relationship matrix A with each matrix mat in diagonal and 0 otherwise
A<-NULL
A<-adiag(mat[[1]],mat[[2]],mat[[3]],mat[[4]],mat[[5]],mat[[6]],mat[[7]],mat[[8]],mat[[9]],mat[[10]],
         mat[[11]],mat[[12]],mat[[13]],mat[[14]],mat[[15]],mat[[16]],mat[[17]],mat[[18]],mat[[19]],mat[[20]],
         mat[[21]],mat[[22]],mat[[23]],mat[[24]],mat[[25]],mat[[26]],mat[[27]],mat[[28]],mat[[29]],mat[[30]],
         mat[[31]],mat[[32]],mat[[33]],mat[[34]])


## Save data
data[[d]]<-list("N"=N,"Y"=Y,"X"=X,
           #L.theta=as.matrix(chol(A)) ## Using decomposition of Cholesky
           'invA'=as.matrix(solve(A)) # inverse of matrix A
)

} #end of loop d


#################
#### ANALYSIS ###
#################
model.fit<-list(); fit.mcmc<-list()
for (d in 1:nSim){ # loop for each simulated datasets in MIXTURE design only!!!

  # ### JAGS ####
  require(R2jags) 
  ### Start of the run ###
  cat( "Sampling MCMC chain using R2JAGS...\n" )
  model.fit[[d]]<- jags(
    model.file = "code/models/LETM.txt",
    data[[d]],inits,parameters,
    n.chains = nChains, n.iter = nPerChain*nChains, n.burnin = burnInSteps,n.thin=thinSteps,
    progress.bar = "text")
  fit.mcmc[[d]] <- as.mcmc(model.fit[[d]])
  
} # end loop d

save(nFather,nMother,elapsed.time,data,inits,model.fit,fit.mcmc,file=paste('results/RESULTS_LETM_SAT_',model.name,'.RData',sep=""))


