rm(list=ls(all=TRUE))


setwd('~/Documents/RESEARCH/PROJECT/R/')
source("HDIofMCMC.R")
source("MCMCsummary.R")
source("functionsgraphmac.R")

library(rjags)

dest<-c('~/Documents/RESEARCH/PROJECT/LETM/BUZATTO/Data for Mat/EARWIGS/')
#dest<-c('~/Dropbox/research2/project/letm0/buzatto/')
setwd(dest)

plotCI <- function(par,n,mcmc,pch){
  res50=as.numeric(mcmc[par,"median"])
  res2.5=as.numeric(mcmc[par,"2.5%"])
  res97.5=as.numeric(mcmc[par,"97.5%"])
  #res25=mcmc[par,"25%"]
  #res75=mcmc[par,"75%"]
  #xrange=c(1,n)
  #yrange=c(min(res2.5),max(res97.5))
  points(n,res50,pch=pch,col=1)
  segments(n,res2.5,n,res97.5,col=1)
  #segments(n,res25,n,res75,lwd=2,col=1)
  #mtext(par, side = 1, line=1, at = n)
  #mtext(c("fullsibs","halfsibs"), side = 1, line=1, at = 1:n,col=1:n)
  #mtext(c("bothdiet","diffdiet"), side = 1, line=1, at = 1:n,col=1:n)
}

# Adds a highest density ellipse to an existing plot
# xy: A matrix or data frame with two columns.
#     If you have to variables just cbind(x, y) them.
# coverage: The percentage of points the ellipse should cover
# border: The color of the border of the ellipse, NA = no border
# fill: The filling color of the ellipse, NA = no fill
# ... : Passed on to the polygon() function
add_hd_ellipse <- function(xy, coverage, border = "blue", fill = NA, ...) {
  library(MASS)
  library(cluster)
  fit <- cov.mve(xy, quantile.used = round(nrow(xy) * coverage))
  points_in_ellipse <- xy[fit$best, ]
  ellipse_boundary <- predict(ellipsoidhull(points_in_ellipse))
  polygon(ellipse_boundary, border=border, col = fill, ...)
}


newcol=gray.colors(3, start = 0.3, end = 0.9, gamma = 2.2, alpha = 0.8)




########################
### EMPIRICAL ##########
########################

#####################
codaSamples=sumMCMC=NULL
#setwd("~/Dropbox/research2/project/letm0/buzatto/Data for Mat/EARWIGS/")
#load('Data for Mat/EARWIGS/RESULTS/RESULTS_LETM_PEDIGREE_all.RData')
load("~/Documents/RESEARCH/PROJECT/LETM/BUZATTO/Data for Mat/EARWIGS/RESULTS/RESULTS_LETM_PEDIGREE_all.RData")

population=c('West Wideopen','Brownsman','East Wideopen')#, 'sta')

data.br <- new.env();load("DATA/br.data.RData",envir=data.br)
data.ewo <- new.env();load("DATA/ewo.data.RData",envir=data.ewo)
data.wwo <- new.env();load("DATA/wwo.data.RData",envir=data.wwo)
#data.sta <- new.env();load("DATA/sta.data.RData",envir=data.sta)

data=list()
data[["wwo"]]<-data.wwo$data
data[["br"]]<-data.br$data
data[["ewo"]]<-data.ewo$data

#data[["sta"]]<-data.sta$data

X <- list(data[[1]]$Pronotum,data[[2]]$Pronotum,data[[3]]$Pronotum)#,data[[4]]$Pronotum)

codaSamples<-fit.mcmc
sumMCMC <- MCMCsummary(codaSamples)
par=rownames(sumMCMC)


table<-cbind(
  median=round(as.numeric(sumMCMC[,"median"]),2),
  q2.5=   round(as.numeric(sumMCMC[,"2.5%"]),2),
  q97.5=       round(as.numeric(sumMCMC[,"97.5%"]),2),
  ROPE=sumMCMC[,"ROPE"])

write.csv(table,file=paste("Figures/Table_Emp.csv",sep=""))


## Tables
MCMCtable(sumMCMC)





load("~/Documents/RESEARCH/PROJECT/LETM/BUZATTO/Data for Mat/EARWIGS/RESULTS/RESULTS_LETM_PEDIGREE_thetas.RData")

N=NULL;nsire=ndam=NULL
for (j in 1:length(population)){
  N[j]=dim(data[[j]])[1]
  #nsire[j]=length(unique(data[[j]]$Sire))
  #ndam[j]=length(unique(data[[j]]$Dam))
}

data_all.tmp = do.call(rbind,data)
data_all <- cbind(data_all.tmp,I=rep(1:length(population),time=N))

I=data_all$I

states<-as.matrix(codaSamples.states)
etas.all=thetas.all=NULL
etas.all <- states[,1:(dim(states)[2]/2)]
thetas.all <- states[,((dim(states)[2]/2)+1) : dim(states)[2]]
# 
# pdf(file = paste('Article/Distributions_medians_Islands.pdf',sep=''),width = 10, height = 8)
# layout(matrix(c(1:4), 2, 2, byrow = TRUE))
# par(mar=c(4.1,4.1,4.1,4.1))
# 
# for (j in 1:length(population)){
# etas=thetas=NULL
#   
#   etas <- etas.all[,which(I==j)]
#   thetas <- thetas.all[,which(I==j)]
#   
#   
#   median.etas <- colMeans(etas) #apply(etas,2,function(x) mean(x))#quantile(x,prob=0.5))
#   median.thetas <- colMeans(thetas) #apply(thetas,2,function(x) mean(X))#quantile(x,prob=0.5))
#   
#   prop.thetas.sup=NULL
#   #for (i in 1:dim(etas)[1]){
#     prop.thetas.sup<-mean(median.thetas > max(median.etas))#/dim(etas)[2]
#   #}
#   
#   smp<-sample(1:dim(etas)[1],100,replace=FALSE)
#   
#   #smp.etas=smp.thetas=array(,dim=c(length(x),dim(etas)[2]))
#   #for (i in 1:length(x)){
#   smp.etas <- etas[smp,]  
#   smp.thetas <- thetas[smp,] 
#   #}
#   
# 
#   X=NULL
#   X <- data[[j]]$Pronotum
#   plot(density(X,adjust=2),type='n', main=population[j],xlim=c(1,3.5),ylim=c(0,6),xlab='Status/Cue')
#   
#   
#   ## Density Threshold
#   theta <- density(median.thetas,adjust=2)
#   #par(new=TRUE);hist(median.thetas, prob=TRUE, col="grey",xlim=c(1,3.5),ylim=c(0,6),breaks=30)
#   polygon(theta,col=newcol[1], border="darkgrey")
#   #rug(df$theta, col="red")
#   #for (i in 1:length(smp)){
#   #  lines(density(smp.thetas[i,]),col=newcol[1])
#   #}
#   #lines(theta,col=1,lwd=3)
#   #hist(theta)
#   
#   ## Density Proximate cue
#   eta <- density(median.etas,adjust=2)
#   polygon(eta, col=newcol[3], border="grey")
#   #par(new=TRUE);hist(median.etas, prob=TRUE, col="grey",xlim=c(1,3.5),ylim=c(0,6),breaks=30)
#   #for (i in 1:length(smp)){
#   #  lines(density(smp.etas[i,]),col=newcol[2])
#   #}
#   #lines(eta,col=2,lwd=3)
#   
#   
#   
#   ## Density Environmental cue
#   #hist(X[[j]], main=population[j],xlim=c(2,3.5),ylim=c(0,10),xlab='Status/Cue',probability=TRUE, border="grey",add=TRUE)
#   dx <- density(X,adjust=2);lines(dx, col=1,lty=2,lwd=3)
#   #rug(X[[j]])
#   
#   # ## Cumulative frequency of Y
#   Y=NULL
#   Y<-ifelse(median.etas>median.thetas,1,0)
#   x=as.numeric(median.etas)
#   out=NULL
#   out <- glm(Y ~ x, family = binomial)
#   
#   invlogit<-function(x) {1/(1+exp(-(x)))}
#   x=seq(1,3.5,by=0.01)
#   p=invlogit(out$coef[1] + out$coef[2]*x)
#   
#   par(new=TRUE)
#   plot(x,p,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1,3.5),lty=3,lwd=2)
#   axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
#   mtext("Probability",side=4,line=2,at=.5)
#   # 
#   # 
#   # for (i in 2:length(smp)){
#   # Y=NULL;Y<-ifelse(smp.etas[i,]>smp.thetas[i,],1,0)
#   # x=NULL;x=smp.etas[i,]
#   # out=NULL;out <- glm(Y ~ x, family = binomial)
#   # x=seq(2,4,by=0.01)
#   # p=NULL;p=invlogit(out$coef[1] + out$coef[2]*x)
#   # lines(x,p,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2.5,3.5),lty=3,lwd=2)
#   # }
#   
#   
# #   legend("topleft",c("Threshold","Proximate cue","Observable cue"),
# #          #title = "Distributions",
# #          box.col="white",bty="n",
# #          col=c(newcol[1:2],1),
# #          #pch=c(15,15,15,NA,NA),
# #          lty=c(1,1,2), lwd=c(1,1,2),
# #          merge=FALSE
# #   )
#   
# legend("topleft",c("Threshold","Proximate cue","Observable cue","Cumulative frequency of Y"),
#        #title = "Distributions",
#        box.col="white",bty="n",
#        col=c(newcol[1:2],1,1),
#        pch=c(15,15,NA,NA),
#        lty=c(1,1,2,3), lwd=c(1,1,2,2)
# )
# 
#  text(3, 5, bquote(P[theta > max(eta)] == .(prop.thetas.sup)))
#   
# }
# 
# dev.off()


######
jpeg("Figures/Figure1.jpeg",res=300,width = 2000, height = 4000)
#pdf(file = paste('Figures/Figure1b.pdf',sep=''),width = 8, height = 10)
layout(matrix(c(1:3), 3, 1, byrow = TRUE))
par(mar=c(2,4,3,3)+0.1,
    oma = c(4,2,1,1) + 0.1)

for (j in 1:length(population)){
  etas=thetas=NULL
  
  etas <- etas.all[,which(I==j)]
  thetas <- thetas.all[,which(I==j)]
  
  
  median.etas <- colMeans(etas) #apply(etas,2,function(x) mean(x))#quantile(x,prob=0.5))
  median.thetas <- colMeans(thetas) #apply(thetas,2,function(x) mean(X))#quantile(x,prob=0.5))
  
  prop.thetas.sup=NULL
  #for (i in 1:dim(etas)[1]){
  prop.thetas.sup<-mean(median.thetas > max(median.etas))#/dim(etas)[2]
  #}
  
  smp<-sample(1:dim(etas)[1],100,replace=FALSE)
  
  #smp.etas=smp.thetas=array(,dim=c(length(x),dim(etas)[2]))
  #for (i in 1:length(x)){
  smp.etas <- etas[smp,]  
  smp.thetas <- thetas[smp,] 
  #}
  
  
  X=NULL
  X <- data[[j]]$Pronotum
  plot(density(X,adjust=2),type='n', main=population[j],xlim=c(1,3.5),ylim=c(0,6),xlab='',ylab='')
  mtext("Density",side=2,line=2,at=3,cex=1.5)
  
  ## Density Threshold
  #theta <- density(median.thetas,adjust=2)
  #par(new=TRUE);hist(median.thetas, prob=TRUE, col="grey",xlim=c(1,3.5),ylim=c(0,6),breaks=30)
  #polygon(theta,col=newcol[1], border="darkgrey")
  #rug(df$theta, col="red")
  #plot(density(smp.thetas[1,]),col=newcol[1])
  for (i in 1:length(smp)){
    lines(density(smp.thetas[i,]),col=newcol[1])
  }
  #lines(theta,col=1,lwd=3)
  #hist(theta)
  
  ## Density Proximate cue
  #eta <- density(median.etas,adjust=2)
  #polygon(eta, col=newcol[3], border="grey")
  #par(new=TRUE);hist(median.etas, prob=TRUE, col="grey",xlim=c(1,3.5),ylim=c(0,6),breaks=30)
  #plot(density(smp.etas[1,]),col=newcol[2])
  for (i in 1:length(smp)){
    lines(density(smp.etas[i,]),col=newcol[2])
  }
  #lines(eta,col=2,lwd=3)
  
  
  
  ## Density Environmental cue
  #hist(X[[j]], main=population[j],xlim=c(2,3.5),ylim=c(0,10),xlab='Status/Cue',probability=TRUE, border="grey",add=TRUE)
  #dx <- density(X,adjust=2);lines(dx, col=1,lty=2,lwd=3)
  #rug(X[[j]])
  
  # ## Cumulative frequency of Y
  #   Y=NULL
  #   Y<-ifelse(median.etas>median.thetas,1,0)
  #   x=as.numeric(median.etas)
  #   out=NULL
  #   out <- glm(Y ~ x, family = binomial)
  #   
  #   invlogit<-function(x) {1/(1+exp(-(x)))}
  #   x=seq(1,3.5,by=0.01)
  #   p=invlogit(out$coef[1] + out$coef[2]*x)
  #   
  #   par(new=TRUE)
  #   plot(x,p,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1,3.5),lty=3,lwd=2)
  #   axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
  #   mtext("Probability",side=4,line=2,at=.5)
  
  
  mu.eta=sd.eta=NULL
  mu.theta=as.numeric(sumMCMC[paste('mu_theta[',j,']',sep=""),'median'])
  var.theta=as.numeric(sumMCMC[paste('sigma2_theta[',j,']',sep=""),'median'])
  mu.eta=as.numeric(sumMCMC[paste('mu_eta[',j,']',sep=""),'median'])
  var.eta=as.numeric(sumMCMC[paste('sigma2_eta[',j,']',sep=""),'median'])
  
  par(new=TRUE)
  x.pred <- seq(1,4,0.01)
  x<-X
  p.pred<-pnorm((x.pred - mu.theta)/sqrt(var.eta[1] + var.theta))
  plot(x.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1,3.5),ylim=c(0,1),lty=1,lwd=1,col="lightgrey")
  #lines(x.pred,p.pred.lower,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
  #lines(x.pred,p.pred.upper,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
  x.obs <-seq(min(x),max(x),0.01)
  p.obs<-pnorm((x.obs - mu.theta)/sqrt(var.eta[1] + var.theta))
  lines(x.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1,3.5),ylim=c(0,1),lty=1,lwd=3)
  
  #   plot(X.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
  #   lines(X.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=3)
  axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
  mtext("Probability",side=4,line=3,at=.5,cex=1.5)
  
  # 
  # 
  # for (i in 2:length(smp)){
  # Y=NULL;Y<-ifelse(smp.etas[i,]>smp.thetas[i,],1,0)
  # x=NULL;x=smp.etas[i,]
  # out=NULL;out <- glm(Y ~ x, family = binomial)
  # x=seq(2,4,by=0.01)
  # p=NULL;p=invlogit(out$coef[1] + out$coef[2]*x)
  # lines(x,p,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2.5,3.5),lty=3,lwd=2)
  # }
  
  
  #   legend("topleft",c("Threshold","Proximate cue","Observable cue"),
  #          #title = "Distributions",
  #          box.col="white",bty="n",
  #          col=c(newcol[1:2],1),
  #          #pch=c(15,15,15,NA,NA),
  #          lty=c(1,1,2), lwd=c(1,1,2),
  #          merge=FALSE
  #   )
  
  legend("topleft",c("Threshold","Proximate cue",
                     #"Observable cue",
                     "Cumulative frequency of Y"),
         #title = "Distributions",
         box.col="white",bty="n",
         col=c(newcol[1:2],1,1),
         pch=c(15,15,NA),
         lty=c(1,1,1), lwd=c(1,1,2)
  )
  
  text(3, 5, bquote(P[theta > max(eta)] == .(prop.thetas.sup)))
  
}
mtext("Status/Cue",side=1,line=3,at=2.25,cex=1.5)

dev.off()






# jpeg("Figures/Figure1.jpeg",res=300,width = 2000, height = 4000)
# #png(file = paste('Figures/Figure1.png',sep=''),res=300,width = 300, height = 700)
# #pdf(file = paste('Figures/Figure1.pdf',sep=''),width = 6, height = 10)
# layout(matrix(c(1:3), 3, 1, byrow = TRUE))
# par(mar=c(2,4,3,2)+0.1,
#     oma = c(4,2,1,1) + 0.1)
# 
# newcol=gray.colors(3, start = 0.3, end = 0.9, gamma = 2.2, alpha = 0.8)
# 
# for (j in 1:3){
#   
#   mu.eta=sd.eta=NULL
#   mu.theta=as.numeric(sumMCMC[paste('mu_theta[',j,']',sep=""),'median'])
#   var.theta=as.numeric(sumMCMC[paste('sigma2_theta[',j,']',sep=""),'median'])
#   mu.eta=as.numeric(sumMCMC[paste('mu_eta[',j,']',sep=""),'median'])
#   var.eta=as.numeric(sumMCMC[paste('sigma2_eta[',j,']',sep=""),'median'])
#   
#   n=10000
#   df <- data.frame(
#     #Group = rep( name.data, each=n), 
#     theta = c(rnorm(n,mean=mu.theta,sd=sqrt(var.theta))),
#     eta = c(rnorm(n,mean=mu.eta,sd=sqrt(var.eta)))
#   )
#   
#   
#   plot(density(X[[j]]),type='n', main=population[j],xlim=c(1.5,3.0),ylim=c(0,10),xlab='',ylab="",cex=1)
#   mtext("Density",side=2,line=2,at=5,cex=1.5)
#   
#   ## Density Threshold
#   theta <- density(df$theta)
#   polygon(theta,col=newcol[1], border="grey")
#   #rug(df$theta, col="red")
#   
#   ## Density Proximate cue
#   eta <- density(df$eta)
#   polygon(eta, col=newcol[2], border="grey")
#   
#   ## Density Environmental cue
#   #hist(X[[j]], main=population[j],xlim=c(2,3.5),ylim=c(0,10),xlab='Status/Cue',probability=TRUE, border="grey",add=TRUE)
#   #d <- density(X[[j]]);lines(d, col=1,lty=2,lwd=2)
#   #rug(X[[j]])
#   
#   ## Cumulative frequency of Y
#   #   Y<-ifelse(df$eta>df$theta,1,0)
#   #   x=df$eta
#   #   out <- glm(Y ~ x, family = binomial)
#   #   
#   #   invlogit<-function(x) {1/(1+exp(-(x)))}
#   #   x=seq(1,4,by=0.01)
#   #   p=invlogit(out$coef[1] + out$coef[2]*x)
#   #   
#   #   par(new=TRUE)
#   #   plot(x,p,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1.5,3),lty=3,lwd=2)
#   #   axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
#   #   mtext("Probability",side=4,line=2,at=.5)
#   
#   ## Cumulative frequency of Y
#   
#   #   Y<-ifelse(median.etas>median.thetas,1,0)
#   #   x=median.etas
#   #   out <- glm(Y ~ x, family = binomial)
#   #   
#   #   invlogit<-function(x) {1/(1+exp(-(x)))}
#   #   x=seq(1,4,by=0.001)
#   #   p=invlogit(out$coef[1] + out$coef[2]*x)
#   #   
#   par(new=TRUE)
#   x.pred <- seq(2,4,0.01)
#   x<-X[[j]]
#   p.pred<-pnorm((x.pred - mu.theta)/sqrt(var.eta[1] + var.theta))
#   plot(x.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1.5,3.0),ylim=c(0,1),lty=1,lwd=1,col="lightgrey")
#   #lines(x.pred,p.pred.lower,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
#   #lines(x.pred,p.pred.upper,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
#   x.obs <-seq(min(x),max(x),0.01)
#   p.obs<-pnorm((x.obs - mu.theta)/sqrt(var.eta[1] + var.theta))
#   lines(x.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(1.5,3.0),ylim=c(0,1),lty=1,lwd=3)
#   
#   #   plot(X.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
#   #   lines(X.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=3)
#   axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
#   mtext("Probability",side=4,line=2,at=.5,cex=1.5)
#   
#   
#   legend("topleft",c("Threshold","Proximate cue",
#                      #"Observable cue",
#                      "Cumulative frequency of Y"),
#          #title = "Distributions",
#          box.col="white",bty="n",
#          col=c(newcol[1:2],1,1),
#          pch=c(15,15,NA),
#          lty=c(NA,NA,1), lwd=c(NA,NA,3),
#          merge=FALSE
#   )
#   
# }
# mtext("Status/Cue",side=1,line=2,at=2.25,cex=1.5)
# dev.off()



#########
#pdf(file = paste('Figures/Figure2.pdf',sep=''),width = 10, height = 8)
jpeg(file=paste('Figures/Figure2.jpeg',sep=''), res=300, width=2000, height=4000)
#layout(matrix(c(1:6), 3, 2, byrow = TRUE))
layout( matrix(c(1,2,2,3,4,4,5,6,6),3,3,T) );layout.show(6)
par(mar=c(4.1,4.1,4.1,4.1))
cor<-array(,dim=c(dim(etas.all)[1],length(population)))
for (j in 1:(length(population))){
  etas=thetas=NULL
  
  etas <- etas.all[,which(I==j)]
  thetas <- thetas.all[,which(I==j)]
  
  for (i in 1:dim(etas.all)[1]){
    cor[i,j]<-cor.test(etas[i,],thetas[i,])$estimate
  }
  
  smp<-sample(1:dim(etas)[1],1,replace=FALSE)
  
  smp.etas <- smp.thetas <- NULL
  smp.etas <- etas[smp,]  
  smp.thetas <- thetas[smp,] 
  smp.cor <- cor[smp,]
  
  
  
  #### Plot 1
  niceplot(nmY=expression(paste("Correlation (",rho,")",sep="")),nmX="", limX=c(0,2), limY=c(-.2,.2),nmlabX=NA,labX=NA)  
  pointbar(1,
           quantile(cor[,j],prob=c(0.5)),
           quantile(cor[,j],prob=c(0.025)),
           quantile(cor[,j],prob=c(0.975)), 
           pch=20,
           colp=1,bgp=1,colb=1,lgt=0.03)
  abline(h=0,lty=2)
  title(main=population[j],cex.main=1.5,col.main="black") 
  
  
  ###### Plot 2
  df <- data.frame(smp.etas,smp.thetas)
  
  ## Use densCols() output to get density at each point
  x <- densCols(smp.etas,smp.thetas, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
                              "#FCFF0070", "#FF940070", "#FF310070"))(256)
  #df$col <- cols[df$dens]
  col<-"#63636370"
  
  niceplot(nmX=expression(paste("Proximate cue (",eta,")",sep="")),nmY=expression(paste("Thresholds (",theta,")",sep="")), 
           limX=setrge(smp.etas), limY=setrge(smp.thetas))  
  title(main=population[j],cex.main=1.5,col.main="black")
  ## Plot it, reordering rows so that densest points are plotted on top
  points(smp.thetas~smp.etas, data=df[order(df$dens),], pch=20, col=col, cex=1.5)
  #,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  
  d<-cbind(smp.etas,smp.thetas)
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636370", lwd=3)
  #legend("topleft", "95% HDI", col = "#63636370", lty = 1, lwd = 3,box.col="white",bty="n")
  
  #with(df, dataEllipse(smp.etas, smp.thetas,  level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label="", center.pch="+",col=newcol[1]))
  
  #library(car)
  #plot(smp.etas,smp.thetas,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  #dataEllipse(smp.etas, smp.thetas, pch=16,levels=c(0.975),col=c(1:2),xlim=c(1.5,2.8),ylim=c(1.5,2.8),xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  
  #require(psych)
  #bvn<-data.frame(V1=smp.etas,V2=smp.thetas)
  #scatter.hist(x=bvn$V1, y=bvn$V2, density=TRUE, ellipse=TRUE,xlim=c(1.5,2.8),ylim=c(1.5,2.8),xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), title=population[j])
  
  #text(round(range(smp.etas),1)[2]-0.2, round(range(smp.thetas),1)[2]-0.05, bquote(rho == .(round(cor.test(smp.etas,smp.thetas)$estimate,2))),cex=2,col=2)
  #legend("topright",c("95% HDI"),col=2,lty=c(1),lwd=c(2),box.col="white",bty="n")
  
  
}


graphics.off()
#dev.off()



#########
#pdf(file = paste('Figures/Figure3.pdf',sep=''),width = 10, height = 8)
jpeg(file=paste('Figures/FigureS1.jpeg',sep=''), res=300, width=2000, height=4000)
#layout(matrix(c(1:6), 3, 2, byrow = TRUE))
layout( matrix(c(1,2,2,3,4,4,5,6,6),3,3,T) );layout.show(6)
par(mar=c(4.1,4.1,4.1,4.1))
cor<-array(,dim=c(dim(etas.all)[1],length(population)))
for (j in 1:(length(population))){
  etas=thetas=NULL
  X <- data_all$Pronotum[which(I==j)]
  etas <- etas.all[,which(I==j)]
  #thetas <- thetas.all[,which(I==j)]
  
  for (i in 1:dim(etas.all)[1]){
    cor[i,j]<-cor.test(etas[i,],X)$estimate
  }
  
  smp<-sample(1:dim(etas)[1],1,replace=FALSE)
  
  smp.etas <- smp.thetas <- NULL
  smp.etas <- etas[smp,]  
  #smp.thetas <- thetas[smp,] 
  smp.cor <- cor[smp,]
  
  
  
  #### Plot 1
  niceplot(nmY=expression(paste("Correlation (",rho,")",sep="")),nmX="", limX=c(0,2), limY=c(0,1),nmlabX=NA,labX=NA)  
  pointbar(1,
           quantile(cor[,j],prob=c(0.5)),
           quantile(cor[,j],prob=c(0.025)),
           quantile(cor[,j],prob=c(0.975)), 
           pch=20,
           colp=1,bgp=1,colb=1,lgt=0.03)
  abline(h=0,lty=2)
  title(main=population[j],cex.main=1.5,col.main="black") 
  
  
  ###### Plot 2
  df <- data.frame(smp.etas,X)
  
  ## Use densCols() output to get density at each point
  x <- densCols(smp.etas,X, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
                              "#FCFF0070", "#FF940070", "#FF310070"))(256)
  #df$col <- cols[df$dens]
  col<-"#63636370"
  
  niceplot(nmX=expression(paste("Proximate cue (",eta,")",sep="")),nmY="Observable cue (X)", 
           limX=setrge(smp.etas), limY=setrge(X))  
  title(main=population[j],cex.main=1.5,col.main="black")
  ## Plot it, reordering rows so that densest points are plotted on top
  points(X~smp.etas, data=df[order(df$dens),], pch=20, col=col, cex=1.5)
  #,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  
  
  #library(car)
  #plot(smp.etas,smp.thetas,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  #dataEllipse(smp.etas, smp.thetas, pch=16,levels=c(0.975),col=c(1:2),xlim=c(1.5,2.8),ylim=c(1.5,2.8),xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  #with(df, dataEllipse(smp.etas, X,  level=0.95, fill=TRUE, fill.alpha=0.2, plot.points=FALSE, add=TRUE,  ellipse.label="", center.pch="+",col=c("#63636370","#63636305")))
  d<-cbind(smp.etas,X)
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636370", lwd=3)
  #legend("topleft", "95% HDI", col = "#63636370", lty = 1, lwd = 3,box.col="white",bty="n")
  
}


graphics.off()
#dev.off()

# 
# ## Diffrences
# 
# delta=seq(-0.1,0.1,length=4)
# yrange=c(0.5,6.5)
# xrange=c(-1,1)
# I=paste("Diff_h2[",1:6,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences")
#      abline(v=0,lty=2)
#        for (i in 1:length(I)){plotCIh(I[i],i,sumMCMC,pch=1)}
# 
# mtext(c(
#          expression(h[br]^2 - h[ewo]^2),
#          expression(h[br]^2 - h[wwo]^2),
#          expression(h[br]^2 - h[sta]^2),
#          expression(h[ewo]^2 - h[wwo]^2),
#          expression(h[ewo]^2 - h[sta]^2),
#          expression(h[wwo]^2 - h[sta]^2)),
#        side = 2, line=0, at = 1:6,las=1)
# legend("topright",population,pch=1:4,box.col="white",bty="n",title="Populations")
#    
# 
# png(file = paste('Article/Differences_Islands.png',sep=''),width = 580, height = 580)
# layout(matrix(c(1:4), 2, 2, byrow = TRUE))
# par(mar=c(4.1,4.1,4.1,4.1))
# 
# delta=seq(-0.1,0.1,length=4)
# yrange=c(0.5,6.5)
# xrange=c(-.5,.5)
# I=paste("Diff_mu_theta[",1:6,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences")
#      abline(v=0,lty=2)
#        for (i in 1:length(I)){plotCIh(I[i],i,sumMCMC,pch=1)}
# 
# mtext(c(
#          expression(mu[theta[br]] - mu[theta[ewo]]),
#          expression(mu[theta[br]] - mu[theta[wwo]]),
#          expression(mu[theta[br]] - mu[theta[sta]]),
#          expression(mu[theta[ewo]] - mu[theta[wwo]]),
#          expression(mu[theta[ewo]] - mu[theta[sta]]),
#          expression(mu[theta[wwo]] - mu[theta[sta]])),
#        side = 2, line=0, at = 1:6,las=1)
# 
# 
# 
# delta=seq(-0.1,0.1,length=4)
# yrange=c(0.5,6.5)
# xrange=c(-0.1,0.1)
# I=paste("Diff_mu_eta[",1:6,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences")
#      abline(v=0,lty=2)
#        for (i in 1:length(I)){plotCIh(I[i],i,sumMCMC,pch=1)}
# 
# mtext(c(
#          expression(mu[eta[br]] - mu[eta[ewo]]),
#          expression(mu[eta[br]] - mu[eta[wwo]]),
#          expression(mu[eta[br]] - mu[eta[sta]]),
#          expression(mu[eta[ewo]] - mu[eta[wwo]]),
#          expression(mu[eta[ewo]] - mu[eta[sta]]),
#          expression(mu[eta[wwo]] - mu[eta[sta]])),
#        side = 2, line=0, at = 1:6,las=1)
#    
# 
# 
# 
# delta=seq(-0.1,0.1,length=4)
# yrange=c(0.5,6.5)
# xrange=c(-0.1,0.1)
# I=paste("Diff_sigma2_theta[",1:6,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences")
#      abline(v=0,lty=2)
#        for (i in 1:length(I)){plotCIh(I[i],i,sumMCMC,pch=1)}
# 
# mtext(c(
#          expression(sigma[theta[br]]^2 - sigma[theta[ewo]]^2),
#          expression(sigma[theta[br]]^2 - sigma[theta[wwo]]^2),
#          expression(sigma[theta[br]]^2 - sigma[theta[sta]]^2),
#          expression(sigma[theta[ewo]]^2 - sigma[theta[wwo]]^2),
#          expression(sigma[theta[ewo]]^2 - sigma[theta[sta]]^2),
#          expression(sigma[theta[wwo]]^2 - sigma[theta[sta]]^2)),
#        side = 2, line=0, at = 1:6,las=1)
# 
# 
# 
# delta=seq(-0.1,0.1,length=4)
# yrange=c(0.5,6.5)
# xrange=c(-0.1,0.1)
# I=paste("Diff_sigma2_eta[",1:6,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences")
#      abline(v=0,lty=2)
#        for (i in 1:length(I)){plotCIh(I[i],i,sumMCMC,pch=1)}
# 
# mtext(c(
#          expression(sigma[eta[br]]^2 - sigma[eta[ewo]]^2),
#          expression(sigma[eta[br]]^2 - sigma[eta[wwo]]^2),
#          expression(sigma[eta[br]]^2 - sigma[eta[sta]]^2),
#          expression(sigma[eta[ewo]]^2 - sigma[eta[wwo]]^2),
#          expression(sigma[eta[ewo]]^2 - sigma[eta[sta]]^2),
#          expression(sigma[eta[wwo]]^2 - sigma[eta[sta]]^2)),
#        side = 2, line=0, at = 1:6,las=1)
# 
# 
# dev.off()
# 









# stat.cor=array(,dim=c(4,3))
# pos.cor=NULL
# for (j in 1:length(population)){
#   #densplot(as.mcmc(cor[,j]), main=population[j])
#   stat.cor[j,]<-quantile(as.mcmc(cor[,j]),prob=c(0.025,0.5,0.975))
#   pos.cor[j]<-mean(as.mcmc(cor[,j])>0)
# }
# 
# table<-cbind(stat.cor,pos.cor);colnames(table)<-c("2.5%","50%","97.5%", "p.pos")
# 


## Cues
#library(car)
#png(file = paste('Article/Correlations.png',sep=''),width = 680, height = 480)
#pdf(file = paste('Article/Correlations.pdf',sep=''),width = 8, height = 6)
jpeg("Figures/FigureS1b.jpeg",res=300,width = 3000, height = 3000)

## Map densities to colors
cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
                            "#FCFF0070", "#FF940070", "#FF310070"))(256)
#df$col <- cols[df$dens]
#col<-c("#FF310070","#63636370")

layout(matrix(c(1:4), 2, 2, byrow = TRUE))
for (j in 1:2){ 
  
  X <- data[[j]]$Pronotum
  n =length(X)
  
  mu.eta=sd.eta=mu.X=sd.X=NULL
  
  mu.theta=as.numeric(sumMCMC[paste('mu_theta[',j,']',sep=""),'median'])
  var.theta=as.numeric(sumMCMC[paste('sigma2_theta[',j,']',sep=""),'median'])
  mu.eta=as.numeric(sumMCMC[paste('mu_eta[',j,']',sep=""),'median'])
  var.eta=as.numeric(sumMCMC[paste('sigma2_eta[',j,']',sep=""),'median'])
  
  eta=NULL
  for (i in 1:length(X)) {eta[i] <- rnorm(1,X[i],sqrt(var.eta))}
  
  df <- data.frame(
    x=X, 
    z = c(rnorm(n,mean=mu.theta,sd=sqrt(var.theta))),
    y = eta
  )
  
  
  
  col<-c("#63636370")#,"#FF310070")
  
  plot(df$y,df$x,type="n",ylab="Observable cue (X)",xlab=expression(paste("Proximate cue (",eta,")",sep="")),xlim=c(1,3),ylim=c(1.5,3),main=population[j])
  # plot points first
  with(df, points(y,x, col=col, pch=20))
  # then ellipses over the top
  #with(d, dataEllipse(y, x,  level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label="", center.pch="+",col=1))
  d<-cbind(df$y,df$x)
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636370", lwd=3)
  
  legend("topleft",
         as.expression( bquote(rho == .(round(cor(df$x,df$y),2)))),
         cex=1.5,
         text.col=1,
         box.col="white",bty="n"
  )
  
  
  
  
  plot(df$y,df$z,type="n",ylab=expression(paste("Threshold (",theta,")",sep="")),xlab=expression(paste("Proximate cue (",eta,")",sep="")),xlim=c(1,3),ylim=c(1.5,3),main=population[j])
  # plot points first
  with(df, points(y,z, col=col, pch=20))
  # then ellipses over the top
  #with(d, dataEllipse(y, z,  level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label="", center.pch="+",col=1))
  d<-cbind(df$y,df$z)
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636370", lwd=3)
  
  legend("topleft",
         as.expression( bquote(rho == .(round(cor(df$z,df$y),2)))),
         cex=1.5,
         text.col=1,
         box.col="white",bty="n"
  )
}
dev.off()




# ## 3D plots
# install.packages("rgl", dependencies = TRUE)
# library(rgl)
# d <- data.frame(I=I,x=X,y=eta,z=theta)
# plot3d(d$x, d$y, d$z, xlab = "X", ylab = "Y", zlab = "Z", type = "p", col = col, size = 4)







## Distributions of parameters
# 
jpeg(file=paste('Figures/FigureS3.jpeg',sep=''), res=300, width=2000, height=4000)
layout(matrix(c(1:2), 2, 1, byrow = TRUE))


# mu
pch=c(16,17)
delta=c(-.1,.1)
xrange=c(.5,3.5)
yrange=c(1.9,2.5)
I=paste("mu_theta[",1:3,"]",sep="")
J=paste("mu_eta[",1:3,"]",sep="")
plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for means (",mu[eta]," & ",mu[theta],")",sep="")))
for (i in 1:length(I)){
  plotCI(I[i],i+delta[1],sumMCMC,pch=pch[1])
  plotCI(J[i],i+delta[2],sumMCMC,pch=pch[2])
}
segments(1,yrange[1],length(I),yrange[1]);
mtext(population,side=1,line=0,at=1:3);mtext("Populations",side=1,line=2,at=2);
#mtext(expression(mu[theta]), side = 2, line=2, at = 2.3)
legend("topleft",c(expression(paste("means of threshold ",mu[theta],sep="")),expression(paste("means of cue ",mu[eta],sep=""))),pch=pch,box.col="white",bty="n")#,title="Parameters")


# sigma
delta=c(-.1,.1)
xrange=c(.5,3.5)
yrange=c(0,0.05)
I=paste("sigma2_theta[",1:3,"]",sep="")
J=paste("sigma2_eta[",1:3,"]",sep="")
plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for variances (",sigma[eta]^2," & ",sigma[theta]^2,")",sep="")))
for (i in 1:length(I)){
  plotCI(I[i],i+delta[1],sumMCMC,pch=pch[1])
  plotCI(J[i],i+delta[2],sumMCMC,pch=pch[2])
}
segments(1,yrange[1],length(I),yrange[1]);
mtext(population,side=1,line=0,at=1:3);mtext("Populations",side=1,line=2,at=2);
#mtext(expression(mu[theta]), side = 2, line=2, at = 2.3)
legend("topleft",c(expression(paste("variances of threshold ",sigma[theta]^2,sep="")),expression(paste("variances of cue ",sigma[eta]^2,sep=""))),pch=pch,box.col="white",bty="n")#,title="Parameters")

# 
# # h2
# xrange=c(.5,4.5)
# yrange=c(0,1)
# I=paste("h2[",1:4,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for heritability (",h^2,")",sep="")))
# for (i in 1:length(I)){plotCI(I[i],i,sumMCMC,pch=1)}
# segments(1,yrange[1],length(I),yrange[1]);
# mtext(c('br','ewo','wwo', 'sta'),side=1,line=0,at=1:4);mtext("Populations",side=1,line=2,at=2.5);
# 
# 
# #rho
# xrange=c(.5,4.5)
# yrange=c(0,1)
# I=paste("rho[",1:4,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for correlation (",rho,")",sep="")))
# for (i in 1:length(I)){plotCI(I[i],i,sumMCMC,pch=1)}
# segments(1,yrange[1],length(I),yrange[1]);
# mtext(c('br','ewo','wwo', 'sta'),side=1,line=0,at=1:4);mtext("Populations",side=1,line=2,at=2.5);
# 
# 
dev.off()






#####################
## Treatment ########
#####################

population <- c('West Wideopen','Brownsman')
data<-list()

wwo <- new.env()
load("RESULTS/LETM.wwo.diffdiet.mixture.RData", envir=wwo)
sumMCMC.wwo <- MCMCsummary(wwo$codaSamples)
data[[1]] <-read.csv(file=paste(dest,"DATA/I wwo with diet treatment.csv",sep=""))

bsmn <- new.env()
load("RESULTS/LETM.br.diffdiet.mixture.RData", envir=bsmn)
sumMCMC.bsmn <- MCMCsummary(bsmn$codaSamples)
data[[2]] <-read.csv(file=paste(dest,"DATA/I bsmn with diet treatment.csv",sep=""))

MCMC <-list(as.matrix(wwo$codaSamples),as.matrix(bsmn$codaSamples))
sumMCMC<-list(sumMCMC.wwo,sumMCMC.bsmn)

table<-list()
for (j in 1:length(population)){ 
table[[j]]<-cbind(
  median=round(as.numeric(sumMCMC[[j]][,"median"]),2),
     q2.5=   round(as.numeric(sumMCMC[[j]][,"2.5%"]),2),
       q97.5=       round(as.numeric(sumMCMC[[j]][,"97.5%"]),2),
ROPE=sumMCMC[[j]][,"ROPE"])

write.csv(table[[j]],file=paste("Figures/Table_Exp_",population[j],".csv",sep=""))

}

par=rownames(sumMCMC.bsmn)


# # Tables
# table<-list()
# MCMCtable <- function(sumMCMC){
#   sumMCMC.table <- sumMCMC #round(sumMCMC,2)
#   
#   HDI <- paste("[",sumMCMC.table[,"2.5%"],";",sumMCMC.table[,"97.5%"],"]",sep="")
#   sumMCMC.table[,4]<-HDI;colnames(sumMCMC.table)[4] <- "95% HDI"
#   table<-sumMCMC.table[,-5]
#   
#   library(knitr)
#   mcmctable<-kable(table, format = "markdown",  padding = 0)#caption = "Title of the table")
#   
#   return(list(table,mcmctable))
# }
# 
# for (j in 1:length(population)){ 
# table[[j]] <- MCMCtable(sumMCMC[[j]])
# write.csv(table[[j]][[1]],file=paste("Figures/Table_Exp_",population[j],".csv",sep=""))
# 
# }



# 
# ## 3D plots
# install.packages("rgl", dependencies = TRUE)
# library(rgl)
# d <- data.frame(I=I,x=X,y=eta,z=theta)
# plot3d(d$x, d$y, d$z, xlab = "X", ylab = "Y", zlab = "Z", type = "p", col = col, size = 4)

# ## Differences
# plotCIh <- function(par,n,mcmc,pch){
#   res50=mcmc[par,"median"]
#   res2.5=mcmc[par,"2.5%"]
#   res97.5=mcmc[par,"97.5%"]
#   #res25=mcmc[par,"25%"]
#   #res75=mcmc[par,"75%"]
#   #xrange=c(1,n)
#   #yrange=c(min(res2.5),max(res97.5))
#   points(res50,n,pch=pch,col=1)
#   segments(res2.5,n,res97.5,n,col=1)
#   #segments(n,res25,n,res75,lwd=2,col=1)
#   #mtext(par, side = 1, line=1, at = n)
#   #mtext(c("fullsibs","halfsibs"), side = 1, line=1, at = 1:n,col=1:n)
#   #mtext(c("bothdiet","diffdiet"), side = 1, line=1, at = 1:n,col=1:n)
# }
# 
# 
# png(file = paste('Article/Differences_treatments.png',sep=''),width = 680, height = 580)
# layout(matrix(c(1:2), 2, 1, byrow = TRUE))
# par(mar=c(4.1,4.1,4.1,4.1))
# 
# delta=c(-0.1,0.1)
# yrange=c(0.5,5.5)
# xrange=c(-1,1)
# I=c("Diff_h2","Diff_mu_eta","Diff_sigma2_eta","Diff_mu_X","Diff_sigma_X")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences among treatments")
# abline(v=0,lty=2)
# for (j in 1:length(population)){ 
#   for (i in 1:length(I)){plotCIh(I[i],i+delta[j],sumMCMC[[j]],pch=j)}
# }
# 
# mtext(
#   c(
#     expression(h[high]^2 - h[low]^2),
#     expression(mu[eta[high]] - mu[eta[low]]),
#     expression(sigma[eta[high]]^2 - sigma[eta[low]]^2),
#     expression(mu[X[high]] - mu[X[low]]),
#     expression(sigma[X[high]] - sigma[X[low]])), 
# side = 2, line=0, at = 1:5,las=1)
# legend("topright",population,pch=1:2,box.col="white",bty="n",title="Populations")
# 
# 
# 
# 
# Diff_mu_theta <- MCMC[[1]][,"mu_theta"] - MCMC[[2]][,"mu_theta"]
# Diff_sigma2_theta=list()
# Diff_sigma2_theta[[1]] <- MCMC[[1]][,"sigma2_theta[1]"] - MCMC[[2]][,"sigma2_theta[1]"]
# Diff_sigma2_theta[[2]] <- MCMC[[1]][,"sigma2_theta[2]"] - MCMC[[2]][,"sigma2_theta[2]"]
# Diff_mu_eta=list()
# Diff_mu_eta[[1]] <- MCMC[[1]][,"mu_eta[1]"] - MCMC[[2]][,"mu_eta[1]"]
# Diff_mu_eta[[2]] <- MCMC[[1]][,"mu_eta[2]"] - MCMC[[2]][,"mu_eta[2]"]
# Diff_sigma2_eta=list()
# Diff_sigma2_eta[[1]] <- MCMC[[1]][,"sigma2_eta[1]"] - MCMC[[2]][,"sigma2_eta[1]"]
# Diff_sigma2_eta[[2]] <- MCMC[[1]][,"sigma2_eta[2]"] - MCMC[[2]][,"sigma2_eta[2]"]
# 
# I <- list(Diff_mu_theta,Diff_sigma2_theta[[1]],Diff_mu_eta[[1]],Diff_mu_eta[[2]],Diff_sigma2_eta[[1]],Diff_sigma2_eta[[2]])
# J <- c(1,1,2,3,2,3)
# 
# delta=c(-0.1,0.1)
# yrange=c(0.5,6.5)
# xrange=c(-.5,.5)
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',yaxt='n',main="Posterior distributions for differences among populations")
# abline(v=0,lty=2)
# for (i in 1:length(I)){
# points(quantile(I[[i]],prob=0.5),i,pch=J[i])
# segments(quantile(I[[i]],prob=0.025),i,quantile(I[[i]],prob=0.975),i)
# }
# mtext(c(
#   expression(mu[theta[wwo]] - mu[theta[bsmn]]),
#   expression(sigma[theta[wwo]]^2 - sigma[theta[bsmn]]^2),
#   expression(mu[eta[wwo]] - mu[eta[bsmn]]),
#   expression(mu[eta[wwo]] - mu[eta[bsmn]]),
#   expression(sigma[eta[wwo]]^2 - sigma[eta[bsmn]]^2),
#   expression(sigma[eta[wwo]]^2 - sigma[eta[bsmn]]^2)),
#   side = 2, line=0, at = 1:6,las=1)
# legend("topright",c("high","low"),pch=2:3,box.col="white",bty="n",title="Treatments")
# 
# dev.off()
# 
# 
# 
# ## Distributions of parameters
# delta=c(-0.1,0.1)
# png(file = paste('Article/Treatment.png',sep=''),width = 480, height = 580)
# layout(matrix(c(1:3), 3, 1, byrow = FALSE))
# 
# #for (j in 1:length(population)){ 
# 
# # mu
# xrange=c(0.5,3.5)
# yrange=c(2.5,3.5)
# I=c("mu_eta[1]","mu_eta[2]","mu_theta")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for means (",mu[eta]," & ",mu[theta],")",sep="")))#,yaxt='n')
# for (j in 1:length(population)){ 
# for (i in 1:length(I)){plotCI(I[i],i+delta[j],sumMCMC[[j]],pch=j)}
# }
# segments(1+delta[1],yrange[1],2+delta[2],yrange[1]);
# mtext(c("High","Low",""), side = 1, line=0.1, at = 1:length(I), cex=0.75)
# mtext(c(expression(mu[eta]),expression(mu[theta])), side = 1, line=2, at = c(1.5,3))
# legend("topleft",population,pch=1:2,box.col="white",bty="n",title="Populations")
#        
# # sigma
# xrange=c(0.5,3.5)
# yrange=c(0,.1)
# I=c("sigma2_eta[1]","sigma2_eta[2]","sigma2_theta[1]")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for variances (",sigma[eta]^2," & ",sigma[theta]^2,")",sep="")))#,yaxt='n')
# for (j in 1:length(population)){ 
#   for (i in 1:length(I)){plotCI(I[i],i+delta[j],sumMCMC[[j]],pch=j)}
# }
# segments(1+delta[1],yrange[1],2+delta[2],yrange[1]);
# mtext(c("High","Low",""), side = 1, line=0.1, at = 1:length(I), cex=0.75)
# mtext(c(expression(sigma[eta]^2),expression(sigma[theta]^2)), side = 1, line=2, at = c(1.5,3))
# legend("topleft",population,pch=1:2,box.col="white",bty="n",title="Populations")
# 
# # h2 & rho
# xrange=c(0.5,4.5)
# yrange=c(0,1)
# I=c("h2[1]","h2[2]","rho[1]", "rho[2]")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for ",h^2," & ",rho,sep="")))
# for (j in 1:length(population)){ 
#   for (i in 1:length(I)){plotCI(I[i],i+delta[j],sumMCMC[[j]],pch=j)}
# }
# segments(1+delta[1],yrange[1],2+delta[2],yrange[1]);
# segments(3+delta[1],yrange[1],4+delta[2],yrange[1]);
# mtext(c("High","Low","High","Low"), side = 1, line=0.1, at = 1:length(I), cex=0.75)
# mtext(c(expression(h^2),expression(rho)), side = 1, line=2, at = c(1.5,3.5))
# legend("topleft",population,pch=1:2,box.col="white",bty="n",title="Populations")
# 
# #}
# 
# dev.off()
# 



# ## ETM
# #png(file = paste('Figures/Figure2.png',sep=''),width = 450, height = 600)
# #pdf(file = paste('Figures/Figure2.pdf',sep=''),width = 8, height = 10)
# jpeg("Figures/Figure2.jpeg",res=300,width = 2000, height = 3000)
# 
# layout(matrix(c(1:2), 2, 1, byrow = TRUE))
# par(mar=c(2,4,3,2)+0.1,
#     oma = c(4,2,1,1) + 0.1)
#    
# 
# 
# for (j in 1:length(population)){
#   
#   mu.eta.lower=mu.eta.upper=mu.eta=var.eta=mu.X=var.X=NULL
#   var.eta.lower=var.eta.upper=NULL
#   
# mu.theta=as.numeric(sumMCMC[[j]]['mu_theta','median'])
# mu.theta.lower=as.numeric(sumMCMC[[j]]['mu_theta','2.5%'])
# mu.theta.upper=as.numeric(sumMCMC[[j]]['mu_theta','97.5%'])
# 
# var.theta=as.numeric(sumMCMC[[j]]['sigma2_theta[1]','median'])
# var.theta.lower=as.numeric(sumMCMC[[j]]['sigma2_theta[1]','2.5%'])
# var.theta.upper=as.numeric(sumMCMC[[j]]['sigma2_theta[1]','97.5%'])
# 
# mu.eta[1]=as.numeric(sumMCMC[[j]]['mu_eta[1]','median'])
# mu.eta.lower[1]=as.numeric(sumMCMC[[j]]['mu_eta[1]','2.5%'])
# mu.eta.upper[1]=as.numeric(sumMCMC[[j]]['mu_eta[1]','97.5%'])
# 
# var.eta[1]=as.numeric(sumMCMC[[j]]['sigma2_eta[1]','median'])
# var.eta.lower[1]=as.numeric(sumMCMC[[j]]['sigma2_eta[1]','2.5%'])
# var.eta.upper[1]=as.numeric(sumMCMC[[j]]['sigma2_eta[1]','97.5%'])
# 
# mu.eta[2]=as.numeric(sumMCMC[[j]]['mu_eta[2]','median'])
# var.eta[2]=as.numeric(sumMCMC[[j]]['sigma2_eta[2]','median'])
# 
# mu.eta[3]=as.numeric(sumMCMC[[j]]['mu_eta[3]','median'])
# var.eta[3]=as.numeric(sumMCMC[[j]]['sd_eta[3]','median'])^2
# 
# 
# 
# 
# mu.X[1]=as.numeric(sumMCMC[[j]]['mu[1]','median'])
# var.X[1]=as.numeric(sumMCMC[[j]]['sigma[1]','median'])^2
# mu.X[2]=as.numeric(sumMCMC[[j]]['mu[2]','median'])
# var.X[2]=as.numeric(sumMCMC[[j]]['sigma[2]','median'])^2
# 
# name.data=c('high','low')
# n=10000
# df <- data.frame(Group = rep( name.data, each=n), 
#                  X = c(rnorm(n,mean=mu.X[1],sd=sqrt(var.X[1])),
#                            rnorm(n,mean=mu.X[2],sd=sqrt(var.X[2]))),
#                  theta = c(rnorm(n,mean=mu.theta,sd=sqrt(var.theta)),
#                            rnorm(n,mean=mu.theta,sd=sqrt(var.theta))),
#                  eta = c(rnorm(n,mean=mu.eta[1],sd=sqrt(var.eta[1])),
#                          rnorm(n,mean=mu.eta[2],sd=sqrt(var.eta[2])))
# )
# 
# 
# 
# X <- data[[j]]$cue
# 
# 
# plot(density(df$theta),type='n', main=population[j],xlim=c(2,4),ylim=c(0,12),xlab='',ylab='')
# mtext("Density",side=2,line=2,at=6,cex=1.5)
# 
# ## Density Threshold
# theta <- density(df$theta)
# polygon(theta,col=newcol[1], border="grey")
# #rug(df$theta, col="red")
# 
# ## Density Proximate cue
# etalow <- density(df$eta[df$Group=="low"])
# polygon(etalow, col=newcol[2], border="grey")
# 
# etahigh <- density(df$eta[df$Group=="high"])
# polygon(etahigh, col=newcol[3], border="grey")
# 
# #etatot <-density(rnorm(n,mean=mu.eta[3],sd=sd.eta[3]))
# #etatot <- density(df$eta)
# #polygon(etatot, col=NA, border="grey",lty=2,lwd=2)
# 
# 
# ## Density Environmental cue
# #hist(X[[j]], main=population[j],xlim=c(2,3.5),ylim=c(0,10),xlab='Status/Cue',probability=TRUE, border="grey",add=TRUE)
# #d <- density(X);lines(d, col=1,lty=2,lwd=2)
# #rug(X[[j]])
# 
# ## Cumulative frequency of Y
# 
# # p.pred.lower<-pnorm((x.pred - mu.theta.lower)/sqrt(var.eta.lower[1] + var.theta.lower))
# # p.pred.upper<-pnorm((x.pred - mu.theta.upper)/sqrt(var.eta.upper[1] + var.theta.upper))
# 
# 
# 
# 
# #   Y<-ifelse(median.etas>median.thetas,1,0)
# #   x=median.etas
# #   out <- glm(Y ~ x, family = binomial)
# #   
# #   invlogit<-function(x) {1/(1+exp(-(x)))}
# #   x=seq(1,4,by=0.001)
# #   p=invlogit(out$coef[1] + out$coef[2]*x)
# #   
# par(new=TRUE)
# 
# # High
# x.pred <- seq(2,4,0.01)
# x<-data[[j]]$cue[data[[j]]$diet=="high"]
# p.pred<-pnorm((x.pred - mu.theta)/sqrt(var.eta[1] + var.theta))
# plot(x.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
# #lines(x.pred,p.pred.lower,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
# #lines(x.pred,p.pred.upper,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
# x.obs <-seq(min(x),max(x),0.01)
# p.obs<-pnorm((x.obs - mu.theta)/sqrt(var.eta[1] + var.theta))
# lines(x.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=3)
# 
# 
# # p.all<-pnorm((x.pred - mu.theta)/sqrt(var.eta[3] + var.theta))
# # lines(x.pred,p.all,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
# # x.all <- seq(min(X),max(X),0.01)
# # p.all<-pnorm((x.all - mu.theta)/sqrt(var.eta[3] + var.theta))
# # lines(x.all,p.all,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=3)
# 
# # Low
# # x<-data[[j]]$cue[data[[j]]$diet=="low"]
# # x.obs <-seq(min(x),max(x),0.01)
# # p.pred<-pnorm((x.pred - mu.theta)/sqrt(var.eta[2]^2 + var.theta^2))
# # p.obs<-pnorm((x.obs - mu.theta)/sqrt(var.eta[2]^2 + var.theta^2))
# # lines(x.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
# # lines(x.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=3)
# 
# axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
# mtext("Probability",side=4,line=2,at=.5,cex=1.5)
# 
# 
# 
# legend("topleft",c("Threshold","Proximate cue (treatment 'Low')","Proximate cue (treatment 'High')",
#                    #"Observable cue",
#                    "Cumulative frequency of Y"),
#        #title = "Distributions",
#        box.col="white",bty="n",
#        col=c(newcol[1:3],1,1),
#        pch=c(15,15,15,NA),
#        lty=c(NA,NA,NA,1), lwd=c(NA,NA,NA,2),
#        merge=FALSE
#        )
# 
# }
# mtext("Status/Cue",side=1,line=3,at=3,cex=1.5)
# 
# dev.off()
# 





# Delta_mu_eta=Delta_sigma2_eta=NULL
# Delta_mu_eta[1] <- mean(mcmc[[1]][,"mu_eta[1]"] - mcmc[[2]][,"mu_eta[1]"] > 0)
# Delta_mu_eta[2] <- mean(mcmc[[1]][,"mu_eta[2]"] - mcmc[[2]][,"mu_eta[2]"] > 0)
# Delta_sigma2_eta[1] <- mean(mcmc[[1]][,"sigma2_eta[1]"] - mcmc[[2]][,"sigma2_eta[1]"] > 0)
# Delta_sigma2_eta[2] <- mean(mcmc[[1]][,"sigma2_eta[2]"] - mcmc[[2]][,"sigma2_eta[2]"] > 0)
# Delta_mu_theta <- mean(mcmc[[1]][,"mu_theta"] - mcmc[[2]][,"mu_theta"] > 0)
# Delta_sigma2_theta <- mean(mcmc[[1]][,"sigma2_theta[1]"] - mcmc[[2]][,"sigma2_theta[1]"] > 0)
# 




#newcol=gray.colors(6, start = 0.3, end = 0.9, gamma = 2.2, alpha = 0.8)

## Comparison between Pop
mcmc=list()
mcmc[[1]] = as.matrix(wwo$codaSamples)
mcmc[[2]] = as.matrix(bsmn$codaSamples)
par=colnames(mcmc[[1]])

mcmc.states<-list(wwo$codaSamples.states,bsmn$codaSamples.states)

## Distributions of parameters
# 
jpeg(file=paste('Figures/FigureS5.jpeg',sep=''), res=300, width=3000, height=4000)
layout( matrix(c(1,2,3,3,4,5,6,6),2,4,byrow=TRUE) );layout.show(4)
par(mar=c(4.1,4.1,4.1,4.1))
par(mar=c(2,4,3,2)+0.1,
    oma = c(4,2,1,1) + 0.1)

treatment <-c("High", "Low")


for (p in 1:2){

  # mu
  pch=c(16,17)
  delta=c(-.2,0,.2)
  xrange=c(.5,3.5)
  yrange=c(2.5,3.5)
  I=paste("mu_theta",sep="")
  J=paste("mu_eta[",1:2,"]",sep="")
  plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n')#,main=expression(paste("Posterior distributions for means (",mu[eta]," & ",mu[theta],")",sep="")))
  
  plotCI(I,1,sumMCMC[[p]],pch=pch[1])
  for (i in 1:length(J)){
    plotCI(J[i],1+i,sumMCMC[[p]],pch=pch[2])
  #  text(1+i,as.numeric(sumMCMC[[p]][J[i],"97.5%"])+0.05,labels=treatment[i],cex=0.5,font=3)
  }
  legend("topleft",c(expression(paste("means of threshold ",mu[theta],sep="")),expression(paste("means of cue ",mu[eta],sep=""))),pch=pch,box.col="white",bty="n")#,title="Parameters")
  segments(1,yrange[1],3,yrange[1]);
  mtext(c("Both",treatment),side=1,line=0,at=1:3,font=3);mtext("Diet treatments",side=1,line=2,at=2);
  #mtext(expression(mu[theta]), side = 2, line=2, at = 2.3)
  text(3.5,0,LETTERS[1])
  
  # sigma
  delta=c(-.2,0,.2)
  xrange=c(.5,3.5)
  yrange=c(0,0.1)
  I=paste("sigma2_theta[1]",sep="")
  J=paste("sigma2_eta[",1:2,"]",sep="")
  plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n')#,main=expression(paste("Posterior distributions for variances (",sigma[eta]^2," & ",sigma[theta]^2,")",sep="")))

    plotCI(I,1,sumMCMC[[p]],pch=pch[1])
    for (i in 1:length(J)){
      plotCI(J[i],1+i,sumMCMC[[p]],pch=pch[2])
     # text(1+i,as.numeric(sumMCMC[[p]][J[i],"97.5%"])+0.005,treatment[i],cex=0.5,font=3)
    }
    legend("topleft",c(expression(paste("variances of threshold ",sigma[theta]^2,sep="")),expression(paste("variances of cue ",sigma[eta]^2,sep=""))),pch=pch,box.col="white",bty="n")#,title="Parameters")
    segments(1,yrange[1],3,yrange[1]);
    mtext(c("Both",treatment),side=1,line=0,at=1:3,font=3);mtext("Diet treatments",side=1,line=2,at=2);
    #mtext(expression(mu[theta]), side = 2, line=2, at = 2.3)
    text(0.11,0,LETTERS[1])
    
    
    
    states=etas=thetas=NULL
    
    states <- as.matrix(mcmc.states[[p]])
    etas <- states[,1:(dim(states)[2]/2)]
    thetas <- states[,((dim(states)[2]/2)+1) : dim(states)[2]]
    
    median.etas <- colMeans(etas) #apply(etas,2,function(x) quantile(x,prob=0.5))
    median.thetas <- colMeans(thetas) #apply(thetas,2,function(x) quantile(x,prob=0.5))
    
    smp<-sample(1:dim(etas)[1],100,replace=FALSE)
    
    #smp.etas=smp.thetas=array(,dim=c(length(x),dim(etas)[2]))
    #for (i in 1:length(x)){
    smp.etas <- etas[smp,]  
    smp.thetas <- thetas[smp,] 
    #}
    
    
    X=NULL;X <- data[[p]]$cue
    d=NULL;d <- data[[p]]$diet
    plot(density(X,adjust=2),type='n', main=population[p],xlim=c(2,4),ylim=c(0,15),xlab='',ylab='',bty="n")
    text(15.5,2,LETTERS[1])
    mtext("Density",side=2,line=2,at=7.5,cex=1)
    
    ## Density Threshold
    #theta <- density(median.thetas,adjust=2)
    #polygon(theta,col=newcol[1], border="grey")
    #rug(df$theta, col="red")
    #lines(theta,col=1,lwd=3)
    for (i in 1:length(smp)){
      lines(density(smp.thetas[i,],adjust=2),col=newcol[1])
    }
    
    
    
    
    ## Density Proximate cue
    #etalow <- density(median.etas[d=="low"],adjust=2)
    #polygon(etalow, col=newcol[5], border="grey")
    #lines(etalow,col=newcol[3],lwd=3)
    for (i in 1:length(smp)){
      lines(density(smp.etas[i,d=="low"],adjust=2),col=newcol[2])
    }
    
    
    ## Density Proximate cue
    #etahigh <- density(median.etas[d=="high"],adjust=2)
    #polygon(etahigh, col=newcol[3], border="grey")
    #lines(etahigh,col=newcol[1],lwd=3)
    for (i in 1:length(smp)){
      lines(density(smp.etas[i,d=="high"],adjust=2),col=newcol[3])
    }
    
    ## Density Environmental cue
    #hist(X[[j]], main=population[j],xlim=c(2,3.5),ylim=c(0,10),xlab='Status/Cue',probability=TRUE, border="grey",add=TRUE)
    #dx <- density(X,adjust=2);lines(dx, col=1,lty=2,lwd=3)
    #rug(X[[j]])
    
    ## Cumulative frequency of Y
    X.pred <- seq(2,4,0.01)
    X.obs <-seq(min(X),max(X),0.01)
    p.pred<-pnorm((X.pred - median(mcmc[[p]][,"mu_theta"]))/sqrt(median(mcmc[[p]][,"sigma2_eta[1]"]) + median(mcmc[[p]][,"sigma2_theta[1]"])))
    p.obs<-pnorm((X.obs - median(mcmc[[p]][,"mu_theta"]))/sqrt(median(mcmc[[p]][,"sigma2_eta[1]"]) + median(mcmc[[p]][,"sigma2_theta[1]"])))
    
    
    #   Y<-ifelse(median.etas>median.thetas,1,0)
    #   x=median.etas
    #   out <- glm(Y ~ x, family = binomial)
    #   
    #   invlogit<-function(x) {1/(1+exp(-(x)))}
    #   x=seq(1,4,by=0.001)
    #   p=invlogit(out$coef[1] + out$coef[2]*x)
    #   
    par(new=TRUE)
    plot(X.pred,p.pred,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=1,col="lightgrey")
    lines(X.obs,p.obs,type='l',bty='n',xlab='',ylab='',xaxt='n',yaxt='n',xlim=c(2,4),lty=1,lwd=3)
    axis(4, at = c(0,.5,1), labels = TRUE, tick = TRUE, line = 0)
    
    mtext("Probability",side=4,line=2,at=0.5,cex=1)
    
    
    
    legend("topleft",
           c("Threshold","Proximate cue (treatment 'Low')","Proximate cue (treatment 'High')",
             #"Observable cue",
             "Cumulative frequency of Y"),
           #title = "Distributions",
           box.col="white",bty="n",
           col=c(newcol[c(1:3)],1,1),
           #pch=c(15,15,15,NA),
           lty=c(1,1,1,1), lwd=c(2,2,2,2)
    )
    
}
mtext("Status/Cue",side=1,line=3,at=3,cex=1)


# 
# segments(1,yrange[1],length(J),yrange[1]);
# mtext(population,side=1,line=0,at=1:2);mtext("Populations",side=1,line=2,at=1.5);
# #mtext(expression(mu[theta]), side = 2, line=2, at = 2.3)

# 
# # h2
# xrange=c(.5,4.5)
# yrange=c(0,1)
# I=paste("h2[",1:4,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for heritability (",h^2,")",sep="")))
# for (i in 1:length(I)){plotCI(I[i],i,sumMCMC,pch=1)}
# segments(1,yrange[1],length(I),yrange[1]);
# mtext(c('br','ewo','wwo', 'sta'),side=1,line=0,at=1:4);mtext("Populations",side=1,line=2,at=2.5);
# 
# 
# #rho
# xrange=c(.5,4.5)
# yrange=c(0,1)
# I=paste("rho[",1:4,"]",sep="")
# plot(xrange,yrange,bty='n',type='n',xlab='',ylab='',xaxt='n',main=expression(paste("Posterior distributions for correlation (",rho,")",sep="")))
# for (i in 1:length(I)){plotCI(I[i],i,sumMCMC,pch=1)}
# segments(1,yrange[1],length(I),yrange[1]);
# mtext(c('br','ewo','wwo', 'sta'),side=1,line=0,at=1:4);mtext("Populations",side=1,line=2,at=2.5);
# 
# 
#dev.off()




#pdf(file = paste('Figures/Figure2b.pdf',sep=''),width = 8, height = 10)
#jpeg("Figures/Figure3.jpeg",res=300,width = 2000, height = 3000)

#layout(matrix(c(1:2), 2, 1, byrow = TRUE))
# par(mar=c(2,4,3,2)+0.1,
#     oma = c(4,2,1,1) + 0.1)

#for (j in 1:length(population)){

#   prop.thetas.sup=NULL
# for (i in 1:length(smp)){
# prop.thetas.sup[i]<-mean(smp.thetas[i,] > max(smp.etas[i,]))/length(smp.etas[i,])
# }
# hist(prop.thetas.sup,probability=TRUE, xlim=c(0,1))
  
#}
#mtext("Status/Cue",side=1,line=3,at=3,cex=1.5)

dev.off()






#########
#pdf(file = paste('Figures/Figure4.pdf',sep=''),width = 10, height = 8)
jpeg(file=paste('Figures/Figure4.jpeg',sep=''), res=300, width=2000, height=4000)
#layout(matrix(c(1:6), 3, 2, byrow = TRUE))
layout( matrix(c(1,2,2,3,4,4),2,3,T) );#layout.show(4)
par(mar=c(4.1,4.1,4.1,4.1))
#cor<-array(,dim=c(dim(etas)[1],length(population)))
for (j in 1:(length(population))){
  states=etas=thetas=NULL
  
  states <- as.matrix(mcmc.states[[j]])
  etas <- states[,1:(dim(states)[2]/2)]
  thetas <- states[,((dim(states)[2]/2)+1) : dim(states)[2]]
  
  median.etas <- colMeans(etas) #apply(etas,2,function(x) quantile(x,prob=0.5))
  median.thetas <- colMeans(thetas) #apply(thetas,2,function(x) quantile(x,prob=0.5))
  
  
  
  X=NULL;X <- data[[j]]$cue
  I=NULL;I <- data[[j]]$diet
  
  cor <-cor.test(etas[1,],thetas[1,])$estimate
  for (i in 2:dim(etas)[1]){
    cor.tmp <-cor.test(etas[i,],thetas[i,])$estimate
    cor <- c(cor,cor.tmp)
  }
  
  smp<-sample(1:dim(etas)[1],1,replace=FALSE)
  
  smp.etas <- smp.thetas <- NULL
  smp.etas <- etas[smp,]  
  smp.thetas <- thetas[smp,] 
  smp.cor <- cor[smp]
  
  
  
  #### Plot 1
  niceplot(nmY=expression(paste("Correlation (",rho,")",sep="")),nmX="", limX=c(0,2), limY=c(-.2,.2),nmlabX=NA,labX=NA)  
  pointbar(1,
           quantile(cor,prob=c(0.5)),
           quantile(cor,prob=c(0.025)),
           quantile(cor,prob=c(0.975)), 
           pch=20,
           colp=1,bgp=1,colb=1,lgt=0.03)
  abline(h=0,lty=2)
  title(main=population[j],cex.main=1.5,col.main="black") 
  
  
  ###### Plot 2
  df <- data.frame(I,smp.etas,smp.thetas)
  
  ## Use densCols() output to get density at each point
  x <- densCols(smp.etas,smp.thetas, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
                              "#FCFF0070", "#FF940070", "#FF310070"))(256)
  #df$col <- cols[df$dens]
  #col<-"#63636370"
  col<-ifelse(df$I=="high","#63636370","#FF310070")
  
  niceplot(nmX=expression(paste("Proximate cue (",eta,")",sep="")),nmY=expression(paste("Thresholds (",theta,")",sep="")), 
           limX=c(2,4), limY=c(2,4))#setrge(smp.thetas))  
  title(main=population[j],cex.main=1.5,col.main="black")
  ## Plot it, reordering rows so that densest points are plotted on top
  #points(smp.thetas~smp.etas, data=df[order(df$dens),], pch=20, col=col, cex=1.5)
  #,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  with(df, points(smp.etas,smp.thetas, col=col, pch=20))
  
  d<-cbind(smp.etas[I=="high"],smp.thetas[I=="high"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636390", lwd=3)
  
  d<-cbind(smp.etas[I=="low"],smp.thetas[I=="low"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636360", lwd=3)
  #legend("topleft", "95% HDI", col = "#63636370", lty = 1, lwd = 3,box.col="white",bty="n")
  
}


graphics.off()
#dev.off()




#########
#pdf(file = paste('Figures/FigureS2.pdf',sep=''),width = 10, height = 8)
jpeg(file=paste('Figures/FigureS2.jpeg',sep=''), res=300, width=2000, height=4000)
#layout(matrix(c(1:6), 3, 2, byrow = TRUE))
layout( matrix(c(1,2,2,3,4,4),2,3,T) );layout.show(4)
par(mar=c(4.1,4.1,4.1,4.1))

for (j in 1:(length(population))){
  
states=etas=thetas=NULL

states <- as.matrix(mcmc.states[[j]])
etas <- states[,1:(dim(states)[2]/2)]
thetas <- states[,((dim(states)[2]/2)+1) : dim(states)[2]]

median.etas <- colMeans(etas) #apply(etas,2,function(x) quantile(x,prob=0.5))
median.thetas <- colMeans(thetas) #apply(thetas,2,function(x) quantile(x,prob=0.5))



X <- data[[j]]$Pronotum
I <- data[[j]]$diet

cor <-cor.test(etas[1,],X)$estimate
for (i in 2:dim(etas)[1]){
  cor.tmp <-cor.test(etas[i,],X)$estimate
  cor <- c(cor,cor.tmp)
}

smp<-sample(1:dim(etas)[1],1,replace=FALSE)

smp.etas <- smp.thetas <- NULL
smp.etas <- etas[smp,]  
smp.thetas <- thetas[smp,] 
smp.cor <- cor[smp]
  
  
  
  #### Plot 1
  niceplot(nmY=expression(paste("Correlation (",rho,")",sep="")),nmX="", limX=c(0,2), limY=c(0,1),nmlabX=NA,labX=NA)  
  pointbar(1,
           quantile(cor,prob=c(0.5)),
           quantile(cor,prob=c(0.025)),
           quantile(cor,prob=c(0.975)), 
           pch=20,
           colp=1,bgp=1,colb=1,lgt=0.03)
  abline(h=0,lty=2)
  title(main=population[j],cex.main=1.5,col.main="black") 
  
  
  ###### Plot 2
  df <- data.frame(I,smp.etas,X)
  
  ## Use densCols() output to get density at each point
  x <- densCols(smp.etas,X, colramp=colorRampPalette(c("black", "white")))
  df$dens <- col2rgb(x)[1,] + 1L
  
  ## Map densities to colors
  cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
                              "#FCFF0070", "#FF940070", "#FF310070"))(256)
  #df$col <- cols[df$dens]
  #col<-"#63636370"
  col<-ifelse(df$I=="high","#63636370","#FF310070")
  
  niceplot(nmX=expression(paste("Proximate cue (",eta,")",sep="")),nmY="Observable cue (X)", 
           limX=setrge(smp.etas), limY=setrge(X))  
  title(main=population[j],cex.main=1.5,col.main="black")
  ## Plot it, reordering rows so that densest points are plotted on top
  #points(X~smp.etas, data=df[order(df$dens),], pch=20, col=col, cex=1.5)
  #,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  # plot points first
  with(df, points(smp.etas,X, col=col, pch=20))
  
  #library(car)
  #plot(smp.etas,smp.thetas,xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  #dataEllipse(smp.etas, smp.thetas, pch=16,levels=c(0.975),col=c(1:2),xlim=c(1.5,2.8),ylim=c(1.5,2.8),xlab=expression(paste("Proximate cue (",eta,")",sep="")),ylab=expression(paste("Thresholds (",theta,")",sep="")), main=population[j])
  #with(df, dataEllipse(smp.etas, X,  level=0.95, fill=TRUE, fill.alpha=0.2, plot.points=FALSE, add=TRUE,  ellipse.label="", center.pch="+",col=c("#63636370","#63636305")))
  d<-cbind(smp.etas[I=="high"],X[I=="high"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636390", lwd=3)
  
  d<-cbind(smp.etas[I=="low"],X[I=="low"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636360", lwd=3)
  #legend("topleft", "95% HDI", col = "#63636370", lty = 1, lwd = 3,box.col="white",bty="n")
  
}


graphics.off()
#dev.off()











# Cues
library(car)
#png(file = paste('Article/Correlations.png',sep=''),width = 680, height = 480)
#pdf(file = paste('Article/Correlations.pdf',sep=''),width = 8, height = 6)
jpeg("Figures/FigureS2b.jpeg",res=300,width = 3000, height = 3000)

## Map densities to colors
cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
                            "#FCFF0070", "#FF940070", "#FF310070"))(256)
#df$col <- cols[df$dens]
#col<-c("#FF310070","#63636370")

layout(matrix(c(1:4), 2, 2, byrow = TRUE))
for (j in 1:length(population)){ 
  
  X <- data[[j]]$cue
  I=NULL;I <- data[[j]]$diet
  
  mu.eta=sd.eta=mu.X=sd.X=NULL
  
  mu.theta=as.numeric(sumMCMC[[j]]['mu_theta','median'])
  sd.theta=sqrt(as.numeric(sumMCMC[[j]]['sigma2_theta[1]','median']))
  theta <- rnorm(length(X),mu.theta,sd.theta)
  
  mu.eta[1]=as.numeric(sumMCMC[[j]]['mu_eta[1]','median'])
  sd.eta[1]=sqrt(as.numeric(sumMCMC[[j]]['sigma2_eta[1]','median']))
  mu.eta[2]=as.numeric(sumMCMC[[j]]['mu_eta[2]','median'])
  sd.eta[2]=sqrt(as.numeric(sumMCMC[[j]]['sigma2_eta[2]','median']))
  
  eta=NULL
  for (i in 1:length(X)) {eta[i] <- rnorm(1,X[i],sd.eta[I[i]])}
  
  df <- data.frame(I=I,x=X,y=eta,z=theta)
  
  col<-ifelse(df$I=="high","#63636370","#FF310070")
  
  plot(df$y,df$x,type="n",ylab="Observable cue (X)",xlab=expression(paste("Proximate cue (",eta,")",sep="")),xlim=c(2,4),ylim=c(2,4),main=population[j])
  # plot points first
  with(df, points(y,x, col=col, pch=20))
  # then ellipses over the top
  #with(d, dataEllipse(x, y, I, level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label=FALSE, center.pch="+"))
  d<-cbind(df$y[df$I=="high"],df$x[df$I=="high"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636390", lwd=3)
  d<-cbind(df$y[df$I=="low"],df$x[df$I=="low"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636360", lwd=3)
  
  legend("topleft",c(
    as.expression( bquote(rho[high] == .(round(cor(df$x[df$I=="high"],df$y[df$I=="high"]),2)))),
    as.expression( bquote(rho[low] == .(round(cor(df$x[df$I=="low"],df$y[df$I=="low"]),2))))),
    text.col=1:2,
    box.col="white",bty="n"
  )
  
  
  
  plot(df$y,df$z,type="n",ylab=expression(paste("Threshold (",theta,")",sep="")),xlab=expression(paste("Proximate cue (",eta,")",sep="")),xlim=c(2,4),ylim=c(2,4),main=population[j])
  # plot points first
  with(df, points(y,z, col=col, pch=20))
  # then ellipses over the top
  #with(d, dataEllipse(z, y, I, level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label=FALSE, center.pch="+"))
  d<-cbind(df$y[df$I=="high"],df$z[df$I=="high"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636390", lwd=3)
  d<-cbind(df$y[df$I=="low"],df$z[df$I=="low"])
  add_hd_ellipse(d, coverage = 0.95, border=NA, fill = "#63636360", lwd=3)
  
  legend("topleft",c(
    as.expression( bquote(rho[high] == .(round(cor(df$z[df$I=="high"],df$y[df$I=="high"]),2)))),
    as.expression( bquote(rho[low] == .(round(cor(df$z[df$I=="low"],df$y[df$I=="low"]),2))))),
    text.col=1:2,
    box.col="white",bty="n"
  )
}
dev.off()







## Cues
#library(car)
#png(file = paste('Article/Correlations.png',sep=''),width = 680, height = 480)
#pdf(file = paste('Article/Correlations.pdf',sep=''),width = 8, height = 6)
# jpeg("Figures/FigureS2.jpeg",res=300,width = 3000, height = 3000)
# 
# ## Map densities to colors
# cols <-  colorRampPalette(c("#00009970", "#00FEFF70", "#45FE4F70", 
#                             "#FCFF0070", "#FF940070", "#FF310070"))(256)
# #df$col <- cols[df$dens]
# #col<-c("#FF310070","#63636370")
# 
# layout(matrix(c(1:4), 2, 2, byrow = TRUE))
# for (j in 1:length(population)){ 
#   
#   X <- data[[j]]$cue
#   I <-data[[j]]$diet
#   
#   mu.eta=sd.eta=mu.X=sd.X=NULL
#   
#   mu.theta=as.numeric(sumMCMC[[j]]['mu_theta','median'])
#   sd.theta=sqrt(as.numeric(sumMCMC[[j]]['sigma2_theta[1]','median']))
#   theta <- rnorm(length(X),mu.theta,sd.theta)
#   
#   mu.eta[1]=as.numeric(sumMCMC[[j]]['mu_eta[1]','median'])
#   sd.eta[1]=sqrt(as.numeric(sumMCMC[[j]]['sigma2_eta[1]','median']))
#   mu.eta[2]=as.numeric(sumMCMC[[j]]['mu_eta[2]','median'])
#   sd.eta[2]=sqrt(as.numeric(sumMCMC[[j]]['sigma2_eta[2]','median']))
#   
#   eta=NULL
#   for (i in 1:length(X)) {eta[i] <- rnorm(1,X[i],sd.eta[I[i]])}
#   
#   d <- data.frame(I=I,x=X,y=eta,z=theta)
#   
#   col<-ifelse(d$I=="high","#63636370","#FF310070")
#   
#   plot(d$x,d$y,type="n",xlab="Observable cue (X)",ylab=expression(paste("Proximate cue (",eta,")",sep="")),xlim=c(2,4),ylim=c(2,4),main=population[j])
#   # plot points first
#   with(d, points(x,y, col=col, pch=20))
#   # then ellipses over the top
#   with(d, dataEllipse(x, y, I, level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label=FALSE, center.pch="+"))
#   
#   legend("topleft",c(
#     as.expression( bquote(rho[high] == .(round(cor(d$x[d$I=="high"],d$y[d$I=="high"]),2)))),
#     as.expression( bquote(rho[low] == .(round(cor(d$x[d$I=="low"],d$y[d$I=="low"]),2))))),
#     text.col=1:2,
#     box.col="white",bty="n"
#   )
#   
#   
# 
# plot(d$z,d$y,type="n",xlab=expression(paste("Threshold (",theta,")",sep="")),ylab=expression(paste("Proximate cue (",eta,")",sep="")),xlim=c(2,4),ylim=c(2,4),main=population[j])
# # plot points first
# with(d, points(z,y, col=col, pch=20))
# # then ellipses over the top
# with(d, dataEllipse(z, y, I, level=0.95, fill=TRUE, fill.alpha=0.1, plot.points=FALSE, add=TRUE,  ellipse.label=FALSE, center.pch="+"))
# 
# legend("topleft",c(
#   as.expression( bquote(rho[high] == .(round(cor(d$z[d$I=="high"],d$y[d$I=="high"]),2)))),
#   as.expression( bquote(rho[low] == .(round(cor(d$z[d$I=="low"],d$y[d$I=="low"]),2))))),
#   text.col=1:2,
#   box.col="white",bty="n"
# )
# }
# dev.off()




