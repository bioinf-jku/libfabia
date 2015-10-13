
nThreads = as.integer(2)
cyc <- as.integer(100)
if (FALSE) {
source("../common.R")
data <- loadDataFromHDF5("mnist.hdf5")
set.seed(42)
X <- data$trainx[, 1:1000]
p = 10
} else {
set.seed(42)
X = matrix(rnorm(20), ncol=4)
p = as.integer(3)
}

dyn.load("libfabia.so")
set.seed(42)
l=as.integer(ncol(X))
n=as.integer(nrow(X))
eps <- as.double(1e-3)
eps1 <- as.double(2.2204460492503131e-15)
init_lapla <- 1.0
init_psi <- 0.2

XX <- sapply(1:n, function(i) sum(X[i, ] * X[i, ])/l);
XX[XX < eps] = eps
nL <- as.integer(0)
lL <- as.integer(0)
bL <- as.integer(0)
alpha <- as.double(0.1)
non_negative <- as.integer(0)
p <- as.integer(p)
spl = as.double(0.0)
spz <- as.double(0.5)
scale <- as.double(0.0)
lap <- as.double(1.0)


"----------- C version ---------"
set.seed(42)
L <- matrix(rnorm(n*p),nrow=n,ncol=p)
Z <- matrix(rnorm(l*p),nrow=p,ncol=l)
lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)
Psi <- init_psi*XX
dummy = NULL
print(L)
#system.time(res1<- .Call("fabia_r", X,Psi,L,lapla,cyc ,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative, nThreads))
system.time(res <- .C("fabia_cm_d_r", as.integer(p), as.integer(n), as.integer(l),
    X=as.double(X), Psi=as.double(Psi), L=as.double(L), Z=as.double(Z),
    lapla=as.double(lapla), cyc, as.double(alpha), as.double(eps), as.double(spl),
    as.double(spz), as.integer(scale), as.double(lap), as.integer(10), nThreads))
#system.double(.C("fabia_cm_f_r", as.integer(p), as.integer(n), as.integer(l),
#    as.single(X), as.single(Psi), as.single(L), as.single(Z),
#    as.single(lapla), cyc, as.single(alpha), as.single(eps), as.single(spl),
#    as.single(spz), as.integer(scale), as.single(lap), as.integer(10), nThreads, DUP=FALSE))

#q()

set.seed(42)
L <- matrix(rnorm(n*p),nrow=n,ncol=p)
Z <- matrix(rnorm(l*p),nrow=p,ncol=l)
lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)
Psi <- init_psi*XX
library(fabia)
system.time(res2 <- .Call("fabic", X,Psi,L,lapla,cyc ,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative,PACKAGE="fabia"))


cat("max diff L:    \t", max(abs(res1$L - res2$L)), "\n")
cat("max diff Z:    \t", max(abs(res1$Z - res2$E_SX_n)), "\n")
cat("max diff Psi:  \t", max(abs(res1$Psi - res2$Psi)), "\n")
cat("max diff lapla:\t", max(abs(res1$lapla - res2$lapla)), "\n")
#cat("\n\n-------------C return values--------------\n");
#print(lapla)



#system.time(res <- .Call("fabic", X,Psi,L,lapla,cyc ,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative,PACKAGE="fabia"))
cat("\n\n-------------R results--------------\n");
X = as.matrix(read.table("X.csv"))
L = as.matrix(read.table("L.csv"))
Z = as.matrix(read.table("Z.csv"))
lapla <-  matrix(1,nrow=l,ncol=p)
Psi <- 0.2*rep(1, n)

XX <- sapply(1:n, function(i) sum(X[i, ] * X[i, ])/l);
kvect <- as.vector(rep(1,p))
nvect <- as.vector(rep(1,n))
nk_one <- matrix(1,n,p)
nk_zero <- matrix(0,n,p)
kk_zero <- matrix(0,p,p)
kk_one <- diag(p)
epsv<-eps1*kvect
epsn<-eps*nvect
E_SX_n <- matrix(0,p,l)
E_SSXX_n <- list()
epsv<-eps*kvect
sum1<- nk_zero
sum2<- eps*kk_one
for (i in 1:cyc){
	LPsi<-diag(1/Psi)%*%L
	LPsiL<-crossprod(L,LPsi)
 	sum1<- nk_zero
	sum2<- eps*kk_one
	for (j in 1:l){
		laj <- lapla[j,]
		tmp <- chol2inv(chol(LPsiL+diag(laj)))
		x_j <- as.vector(X[,j])
		e_sx_n <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
		sum1 <- sum1 +  tcrossprod(x_j,e_sx_n)
		#cat("\nsum1\n");print(sum1)
		#e_sx_n[which(e_sx_n<0)] <- 0
		e_ssxx_n <-  tmp + tcrossprod(e_sx_n,e_sx_n)
		sum2 <- sum2 + e_ssxx_n
		laj <- (epsv+diag(e_ssxx_n))^(-spz)
		laj[which(laj<lap)] <- lap
		lapla[j,] <- laj
		#print(laj)
	}

    sll <- chol2inv(chol(sum2))
    L <- sum1%*%sll
    ddL <- alpha*Psi*sign(L)*abs(nk_one*eps+L)^{-spl}
    L <- L - ddL
    L[which(abs(L)<abs(ddL))] <- 0
	#cat("\nL\n");	print(L)
	Psi <- pmax(epsn, abs(XX - diag(tcrossprod(L,sum1))/l))
	#cat("\nPsi\n");	print(Psi)
	#print(diag(tcrossprod(L,sum1)))
}

#L
for (j in 1:l){
	laj <- lapla[j,]
	tmp <- chol2inv(chol(LPsiL+diag(laj)))
	x_j <- as.vector(X[,j])
	Z[, j] <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
}

################################################################


if (FALSE) {

l=as.integer(4)
n=as.integer(5)
p = as.integer(3)
cyc = as.integer(100)

if (TRUE) {
	set.seed(42)
	X = matrix(rnorm(l*n), ncol=l)
	L <- matrix(rnorm(n*p),nrow=n,ncol=p)
	Z <- matrix(rnorm(l*p),nrow=p,ncol=l)

	write.table(X, "X.csv", col.names=FALSE, row.names=FALSE)
	write.table(L, "L.csv", col.names=FALSE, row.names=FALSE)
	write.table(Z, "Z.csv", col.names=FALSE, row.names=FALSE)
} else {
	X = as.matrix(read.table("X.csv"))
	L = as.matrix(read.table("L.csv"))
	Z = as.matrix(read.table("Z.csv"))
}

lapla <-  matrix(1,nrow=l,ncol=p)
Psi <- 0.2*rep(1, n)
library(fabia)
alpha <- as.double(0.1)
non_negative <- as.integer(0)
p <- as.integer(p)
spl = as.double(0.0)
spz <- as.double(0.5)
scale <- as.double(0.0)
lap <- as.double(1.0)
eps <- as.double(1e-3)
eps1 <- as.double(2.2204460492503131e-15)
nL <- as.integer(0)
lL <- as.integer(0)
bL <- as.integer(0)

#ref <- .Call("fabic", X,Psi,L,lapla,cyc ,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative,PACKAGE="fabia")



"----------- libfabia version ---------"
dyn.load("libfabia.so")

res <- .C("fabia_cm_f_r", p, n, l, 
	X=as.single(X),Psi=as.single(Psi),L=as.single(L), Z=as.single(Z), lapla=as.single(lapla),
	as.integer(cyc) ,as.single(alpha), as.single(eps), as.single(spl), as.single(spz) ,as.integer(scale),
	as.single(lap), as.integer(1), as.integer(1))


res$L <- matrix(res$L, nrow=n)
res$Z <- matrix(res$Z, nrow=p)
res$lapla <-  matrix(res$lapla, nrow=l)
