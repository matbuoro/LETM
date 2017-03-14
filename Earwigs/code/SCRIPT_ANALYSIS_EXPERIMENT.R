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
#load("EXPERIMENTS.RData") # contains data for 2 populations

# Choose between popualtions "wwo" and "br" to analyze
population='br' #'wwo' #'br' 
data <- read.csv(paste("data/DATA_EXPERIMENTS_",population,".csv",sep=""),sep=";")


##------------- FUNCTION --------##
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
  
  return(list(nAsc=nAsc,nMother=nMother,nFather=nFather,pedigree=as.data.frame(ped)))
}


##--------Build datasets----------##
  
  N=dim(data)[1] # nb of individuals / population

  data.tmp <- data
  data.tmp <- data.tmp[order(data.tmp$family),] # sort by families
  nfamily = length(unique(data.tmp$family)) # nb families
  
  # Generating pedigree (mixture of fullsibs and halfsibs)
  pedigree=list()
  for (n in 1:nfamily){
    nDesc=table(data.tmp$family)[n]
    nMother=1 # unique mother / family 
    nFather=rpois(1,2)+1 # # at least 1 father / family 
    nAsc=nMother + nFather
    data.ped=GeneratePedigree(nDesc,nMother, nFather,n)
    pedigree[[n]]=data.ped$ped
  }

  
  df<-do.call("rbind", pedigree) # dataframe
  df$id=1:length(df$ID) # new index
  
  # new index for ascendants:
  for (i in 1:dim(df)[1]){
    df.tmp = subset(df,df$family==df$family[i])
    df$Sire[i]=ifelse( is.na(df$Sire[i])==T,df$Sire[i], df.tmp$id[df$Sire[i]])
    df$Dam[i]=ifelse( is.na(df$Dam[i])==T,df$Dam[i], df.tmp$id[df$Dam[i]])
  }
  
  pedigree <- data.frame(ID=1:dim(df)[1],Dam=df$Dam,Sire=df$Sire,family=df$family)
  pedigree <- as.data.frame(na.omit(pedigree))  
  
  data<-cbind(data.tmp,pedigree)


## combining all datasets in one
#data_all.tmp = do.call(rbind,data)
#data_all <- cbind(data_all.tmp,I=rep(1:length(population),time=N))


## Re-indexing of the Sire and dam ID for the whole dataset
# here, we add 1000 to make sure all sire and dams of each diet treatments have different ID
## To avoid overlapping IDs between Sire and Dams:   
data$ID <- data$ID + 1000
data$Sire <- data$Sire + 1000
data$Dam <- data$Dam + 1000

sire.tmp <- unique(data$Sire);new.sire.id <- NULL
dam.tmp <- unique(data$Dam);new.dam.id <- NULL

for (i in 1:dim(data)[1]){
  for (j in 1:length(sire.tmp)){
    if(data$Sire[i] == sire.tmp[j]) new.sire.id[i] <- j
  }
  
  for (j in 1:length(dam.tmp)){
    if(data$Dam[i] == dam.tmp[j]) new.dam.id[i] <- length(sire.tmp) + j
  } 
}
data$Sire <- new.sire.id
data$Dam <- new.dam.id



#### 
data <- data[order(data$diet),] # sort by diet / diet treatment "high" go first

Y=data$Morph # alternative phenotypes 1/0
X=data$Pronotum # Observable cue
diet=as.numeric(data$diet) # Indicator of diet treatment

data.jags.all <- list(
  ndiet=max(diet),
  diet=diet,
  n=table(diet), # nb of individus / diet treatments
  X=X,
  Y=Y,
  Ndesc=dim(data)[1], 
  Nasc=max(data$Sire,data$Dam),
  Sire=as.vector(data$Sire),
  Dam=as.vector(data$Dam)
)


#---------------------------LETM MODEL----------------------------------#
write("
  
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
      
      ", "code/models/LETM_EXP.md")

#---------------------------ANALYSIS----------------------------------#


## PARAMETERS TO SAVE
parameters<-c(
  "mu_theta","sigma2_T","sigma2_eta","sigma2_theta","h2","mu_eta","sd_eta",
  "Diff_mu_eta","Diff_sigma2_eta"#,"Diff_sigma2_T","Diff_h2"
)

nchains = 2 # Number of chains to run.
nadapt= 1000
nburnin=5000
nstore=20000 # Total number of steps in chains to save.
nthin=1 # Number of steps to "thin" (1=keep every step).

#===========================
burnin = function(iter){
  if (iter<=10000) {nburnin = iter * 0.2} else { nburnin = 2500};
  return(nburnin)
}

## INITS
mu_theta=mean(X)
eta=NULL
# simualting individual proximate cues
for (i in 1:dim(data)[1]){
  eta[i]=ifelse(data.jags.all$Y[i]==0,
                runif(1,mu_theta-0.1,mu_theta-0.01),
                runif(1,mu_theta+0.01,
                      mu_theta+0.1))
}

inits = function(){list(
  mu_theta=mean(X),
  eta=eta,
  #tau_T=runif(1,0,0.1),
  sigma2_T=c(runif(1,0,10),NA),
  h2=rbeta(2,2,2)
)}



cat( "Adapt...\n" )
model.fit <- jags.model( file = "code/models/LETM_EXP.md", 
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

save(codaSamples,file=paste('results/RESULTS_LETM_EXPERIMENTS_',population,'.RData',sep=""))


# # To save individuals thresholds and cues:
# codaSamples.states <- coda.samples(model.fit,
#                             var = c("theta","eta"),
#                             n.iter = 1000*10,thin = 10)
# 
# save(codaSamples.states,file=paste('results/RESULTS_LETM_EXPERIMENTS_',population,'_States.RData',sep=""))


#==========Outputs=================
# summary(codaSamples)
# library(mcmcplots)
# traplot(codaSamples, "h2")
# caterplot(codaSamples, "h2")

