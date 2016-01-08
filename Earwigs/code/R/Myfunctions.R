#will install the x library if you don't have it
check.packages <- function(x) {
  if (!x %in% installed.packages()){
    install.packages(x);
    cat(paste("Package(s)",x,"installed\n",sep=" "))
  } else {
    cat(paste("Package(s)",x,"already installed\n",sep=" "))
  }
  library(x,character.only = TRUE)
  if(require(x,character.only = TRUE)) cat(paste("Package(s)",x,"loaded\n",sep=" "))
};
#check.packages("Matrix")

CondFactor=function(W,L) {(1e5 * W)/(L^3)}

invlogit<-function(x) {1/(1+exp(-(x)))}

myquantile = function(datavec) { return(quantile(datavec,probs=c(0.05, 0.5, 0.95))) }

Rvector = function(vector) {
  out = paste("c(",paste(x, sep="", collapse=","),")",sep="")
  cat("\n",out,"\n","\n")
  return(out)
}

Rmatrix = function(matrix) {
  out = paste("matrix(", Rvector(as.matrix(mat)), "," ,nrow(mat), ",", ncol(mat), ")", sep="")
  cat("\n",out,"\n","\n")
  return(out)
}


### Calcul du % de valeurs positives #######
## region of practical equivalence (ROPE)
ROPE<-function(x,val){
    # Calcul des % des param?tres
    x.sup<-x[x>val]    # valeurs > ? 0
    x.inf<-x[x<val]    # valeurs < ? 0
    # % de valeurs > ? 0
    prop.x.inf<-round((length(x.inf)/length(x))*100,1)
    prop.x.sup<-round((length(x.sup)/length(x))*100,1)
    rope = paste(prop.x.inf,"% < ",val," < ",prop.x.sup,"%",sep="")
    return(rope)
}





# Tirage des traits suivant une loi normale bivari¬ªe
rbivnorm <- function(N,mux,muy, sigmax,sigmay,rho){
  x <- rnorm(N,mux,sigmax)
  y <- rnorm(N,muy+rho*sigmay*(x-mux)/sigmax, sigmay*sqrt(1-rho^2))
  cbind(x,y)
}

## Calculate n.burnin
burnin = function(iter){
  if (iter<10000) {nburnin = iter * 0.2} else { nburnin = 2500};
  return(nburnin)
};


# Create .md, .html, and .pdf files
Rmdconverter = function(filename,type){
  stopifnot(require(knitr))
  stopifnot(require(markdown))
  # Create .md, .html, and .pdf files
  knit(paste(filename,".Rmd",sep=""))
  markdownToHTML(paste(filename,".md",sep=""), paste(filename,".html",sep=""), options=c("use_xhml"))
  if(type == 'pdf') {system(paste("pandoc -s ",filename,".html -o ",filename,".pdf",sep=""))}
  if(type == 'doc') {system(paste("pandoc -s ",filename,".html -o ",filename,".docx",sep="")) }
}
