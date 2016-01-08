# N: number of families
# T: total number of individuals
# mi : number of individuals of a given morph in familiy i
# ni : number of individuals in familiy i
# h2: heritability measured on the underlying scale
# h2_01: heritability measured on the 0,1 scale

library(MASS)
h2 <- 0.5		# Heritability
Vp <- 1			# Phenotypic variance
Va <- h2 * Vp		# Additive genetic variance
Ve <- (1-h2)*Vp		# Environmental variance
SD.a <- sqrt(Va)	# Additive genetic standard deviation
SD.e <- sqrt(Ve)	# Environmental standard deviation
mu <- 0			# Mean genetic value

N=20; nind=20; n=rep(nind,N); T=sum(n)
K<-rep(1:N,each=nind)
#A <- array(0.5,dim=c(nind,nind))
#diag(A)<-1

# FULL-SIBS
app.coef <- 0.5 	# Relatedness (0.5 for parent offspring pairs and full-siblings, 0.25 for half siblings, 0.125 for first cousins etc).

A <- matrix(0,T,T)		# Additive genetic relationship matrix (relatedness coefficient): A
	sibship=rep(1:N, each = nind)	# Batch index
	for (j in 1:length(sibship)){
	A[,j]<-ifelse(sibship[]==sibship[j],app.coef,0)
	}
	diag(A)<-1
	#fix(A)
	
CovA<-A*Va

Theta <- mvrnorm(1, mu=rep(0,T), Sigma=CovA)

X_tmp <- rnorm(T,0,1)				# Observable cue
X<-scale(X_tmp)
repeat {
Eta<-rnorm(T,X,SD.e)		# Proximate cue
if(max(abs(var(Eta)- (var(X)+Ve)))<abs(0.01)) { break }
	} # end repeat

Y <- ifelse(Eta > Theta, 1, 0)

K<-rep(1:N,each=nind)

m=NULL
for (i in 1:N){m[i]<-sum(Y[K==i])}

## FULL-SIB DESIGN
p=(sum(m/n))/N
c=sqrt(-log(4*p*(1-p)))
x=(sign(0.5-p)*(1.238*c*(1+0.0262*c)))
z=(exp(-0.5*(x^2)))/(sqrt(2*pi))
MSaf<-(sum((m^2)/n) - ((sum(m)^2)/T))/(N-1) # mean squares among families
MSap<-(sum(m) - sum((m^2)/n))/(T-N) # mean squares within families
k= (T- sum((n^2)/T))/(N-1)
t<-(MSaf - MSap)/ (MSaf + (k-1)*MSap)

h2e=2*t*((p*(1-p))/z^2)
h2e_01=2*t

SE_h2e= ((p*(1-p))/z^2)*(2*(1-t))*(1+(k-1)*t)*sqrt((2*(T-1))/((k^2)*(T-N)*(N-1)))


## HALF-SIB DESIGN
Di: number of dams of the ith sires
pij: prpoprtion of a given morph in the family of the jth dam mated to the ith sire 
Dt: total number of dams
S: number of sires
MSap=
MSad=
MSas=

Vap=MSap
Vad=(MSad-MSap)/k1
Vas=(MSas - (MSap+k2*Vad))/k3
Vaf=(MSaf-MSap)/n

p=((sum(sum(pij)))*(1/Di))/S

h2=2*t*((p*(1-p))/z^2)



