GeneratePedigree = function(nDesc,nMother, nFather,id.Family){
  
  nAsc=nMother+nFather # Number of Ascendants
  N = nAsc + nDesc # Total number of individuals
  ped=array(NA,dim=c(N,4));colnames(ped)<-c("ID","Sire","Dam","family") # pedigree matrix
  ped[,"ID"]=1:N # index of individuals
  ped[,"family"]= id.Family # index of family
  id.dam=1:nMother
  id.sire=(nMother+1):(nMother+nFather)
  
  for (i in (nAsc+1):dim(ped)[1]){
    ped[i,"Sire"]=ifelse(nFather==1,id.sire,sample(id.sire,1,replace=TRUE)) 
    ped[i,"Dam"]= ifelse(nMother==1,id.dam,sample(id.dam,1,replace=TRUE))
  }
  
  return(list(nAsc=nAsc,ped=as.data.frame(ped)))
}


SimPedigree = function(nDesc,nFamily,type){

## 1. GENERATE DATASET

P=list();  nMother=nFather=NULL
for (j in 1:nf){
  
  ## sampling of number of father and mother for each batches:
  if(type == "mixture"){
    repeat{
      nMother[j]=1
      nFather[j]=rpois(1,2)
      if(nFather[j]>0 ){break}}  
  }
  if(type == "fullsibs"){nMother[j]=1; nFather[j]=1}
  
  ## Simulate pedigree
  simul.ped=GeneratePedigree(nDesc,nMother[j],nFather[j],j)
  P[[j]]=as.data.frame(simul.ped$ped)
}

df<-do.call("rbind", P) # dataframe
df$id=1:length(df$ID) # new index

# new index for ascendants:
for (i in 1:dim(df)[1]){
 df.tmp = subset(df,df$family==df$family[i])
 df$Sire[i]=ifelse( is.na(df$Sire[i])==T,df$Sire[i], df.tmp$id[df$Sire[i]])
 df$Dam[i]=ifelse( is.na(df$Dam[i])==T,df$Dam[i], df.tmp$id[df$Dam[i]])
}

pedigree <- data.frame(ID=1:dim(df)[1],Dam=df$Dam,Sire=df$Sire,family=df$family)

return(list(pedigree=pedigree,nMother=nMother,nFather=nFather))
}

