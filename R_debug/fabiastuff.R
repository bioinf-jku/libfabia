

m=100
n=1000
l=100
alpha_fabia = 0.4
alpha_fabias = 0.6
spl=1.0
spz=1.0
cyc=200
p=13

library(fabia)
perf <- c()

dat <- makeFabiaData(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,of1 = 5,of2 = 10,  sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,sd_z = 1.0, sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

res <- myfabia(dat$X,p=13, alpha=0.01, center=0, norm=0)
bic = extractBic(res)
sampleList = list()
featureList = list()
for (i in 1:bic$np) {
    sampleList[[i]] = bic$numn[i,]$numnp
    featureList[[i]] = bic$numn[i,]$numng
}
s = consensus_score(n, l, sampleList, featureList, dat$ZC, dat$LC)
print(s)


perf <- c(perf, indrfabia[2])


zz = matrix(0, nrow=10, 100)
for (i in seq_along(dat$ZC))
    zz[i, dat$ZC[[i]]] = 1
ll = matrix(0, nrow=10, 1000)
for (i in seq_along(dat$LC))
    ll[i, dat$LC[[i]]] = 1   

write.table(dat$X, file="/home/tom/X.csv", row.names=F, col.names=F)
write.table(zz, file="/home/tom/row_orig.csv", row.names=F, col.names=F)
write.table(ll, file="/home/tom/col_orig.csv", row.names=F, col.names=F)
write.table(res@L, file="/home/tom/Lres.csv", row.names=F, col.names=F)
write.table(res@Z, file="/home/tom/Zres.csv", row.names=F, col.names=F)



# import paper results
import sklearn.metrics
import numpy as np
Xn = np.loadtxt("/home/tom/X.csv").T
ro = np.loadtxt("/home/tom/row_orig.csv")
co = np.loadtxt("/home/tom/col_orig.csv")
model = FabiaBiclustering(13, alpha = 0.01)
model.fit(Xn)
print sklearn.metrics.consensus_score(model.biclusters_, (ro, co))


def to_list(mat):
    mat = mat.astype(np.bool)
    n = mat.shape[0]
    mat = np.split(mat, n)
    for i in range(n): mat[i] = mat[i].ravel()
    return mat

zc = to_list(zc)
lc = to_list(lc)
fz = to_list(zf)
fl = to_list(lf)
print sklearn.metrics.consensus_score((zf, lf), (zc, lc))








import sklearn.metrics
from fabia import *
import cPickle as pickle

# create data
k, n, m = 10, 100, 1000
random_state = 42
if not os.path.exists("originaldata.pkl"):
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), k, (5, 25),
        (10, 210), 3, 0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
    with open("originaldata.pkl", "w") as f:
        pickle.dump((Xn, X, zc, lc), f, -1)
    
    np.savetxt("/home/tom/X.csv", Xn.T)
    Xn = np.loadtxt("/home/tom/X.csv").T
    
with open("originaldata.pkl", "r") as f:
    (Xn, X, zc, lc) = pickle.load(f)
    
# run python version
model = FabiaBiclustering(k, 400, alpha=0.1, random_state=0)
model.fit(Xn)
print sklearn.metrics.consensus_score(model.biclusters_, (zc, lc))

# parameters from paper
model = FabiaBiclustering(13, 400, alpha=0.4, spl=1.0, spz=1.0, random_state=0)
model.fit(Xn)
print sklearn.metrics.consensus_score(model.biclusters_, (zc, lc))


# import R results
model.thresZ = 0.5
k, n, m = 13, 100, 1000
L = np.loadtxt("/home/tom/Lres.csv").T
Z = np.loadtxt("/home/tom/Zres.csv").T
mom = np.sum((L**2).sum(1) * (Z**2).sum(0)) / (k*n*m)
thresL = np.sqrt(mom) / model.thresZ
model.columns_ = [np.abs(L[i, :]) > thresL for i in range(k)]
model.rows_ = []
for i in range(k):
    idx = np.where(np.abs(Z[:, i]) > model.thresZ)[0]
    sp = np.where(Z[:, i] > model.thresZ, Z[:, i], 0).sum()
    sn = -np.where(Z[:, i] < -model.thresZ, Z[:, i], 0).sum()
    if sp > sn:
        model.rows_.append(Z[:, i] > model.thresZ)
    else:
        model.rows_.append(Z[:, i] < -model.thresZ)
zc = np.loadtxt("/home/tom/row_orig.csv")
lc = np.loadtxt("/home/tom/col_orig.csv")
print sklearn.metrics.consensus_score(model.biclusters_, (zc, lc))


mom <- mom + sum(noL[,i]^2)*sum(nZ[i,]^2)

# import paper results
zc = np.loadtxt("/home/tom/zc.csv")
lc = np.loadtxt("/home/tom/lc.csv")

fz = np.loadtxt("/home/tom/fz.csv")
fl = np.loadtxt("/home/tom/fl.csv")

def to_list(mat):
    mat = mat.astype(np.bool)
    n = mat.shape[0]
    mat = np.split(mat, n)
    for i in range(n):
        mat[i] = mat[i].ravel()
    
    return mat

zc = to_list(zc)
lc = to_list(lc)
fz = to_list(fz[0:8, :])
fl = to_list(fl[0:8, :])
print sklearn.metrics.consensus_score((fz, fl), (zc, lc))

'''

load("exp_1.RData")
data = dat
write.table(dat$X, file="X.csv", row.names=F, col.names=F)
load("exp_1_Biclust_fabia.RData")

# true results
zz = matrix(0, nrow=10, 100)
for (i in seq_along(data$ZC))
    zz[i, dat$ZC[[i]]] = 1

ll = matrix(0, nrow=10, 1000)
for (i in seq_along(data$LC))
    ll[i, dat$LC[[i]]] = 1

write.table(fl, file="lc.csv", row.names=F, col.names=F)
write.table(fz, file="zc.csv", row.names=F, col.names=F)

# fabia results
fl = t(rfabia@RowxNumber)
fz = rfabia@NumberxCol
write.table(fl, file="fl.csv", row.names=F, col.names=F)
write.table(fz, file="fz.csv", row.names=F, col.names=F)

# rerun fabia with paper parameters
X = as.matrix(read.table("X.csv"))
res = fabia(X, p=13, cyc=200, alpha=0.4, spl=1, spz=1)
write.table(res@L, file="L.csv", row.names=F, col.names=F)
write.table(res@Z, file="Z.csv", row.names=F, col.names=F)


library(fabia)

X = as.matrix(read.table("X.csv"))
res = fabia(X, 10, cyc=400)
write.table(res@L, file="L.csv", row.names=F, col.names=F)
write.table(res@Z, file="Z.csv", row.names=F, col.names=F)


data <- makeFabiaData(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,of1 = 5,of2 = 10,
                      sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,sd_z = 1.0,
                      sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)
zz = matrix(0, nrow=10, 100)
for (i in seq_along(data$ZC))
    zz[i, data$ZC[[i]]] = 1

ll = matrix(0, nrow=10, 1000)
for (i in seq_along(data$LC))
    ll[i, data$LC[[i]]] = 1    

write.table(data$X, file="X.csv", row.names=F, col.names=F)
write.table(zz, file="zc.csv", row.names=F, col.names=F)
write.table(ll, file="lc.csv", row.names=F, col.names=F)




data <- makeFabiaData(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,of1 = 5,of2 = 10,
                    sd_noise = 3.0,
                    sd_z_noise = 0.2, mean_z = 2.0,sd_z = 1.0,
                    sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)
res <- fabia(data$X, p=10) 
library(colorspace)
col <- diverge_hcl(25, c = 100, l = c(50, 90), power = 1)
image(t(data$Y), col = col, main = "Data (w/o noise)")

write.table(data$Y, file="Y.csv", row.names=F, col.names=F)
write.table(data$X, file="X.csv", row.names=F, col.names=F)


'''




def robust_scale(X, copy=True, axis=0):
    if copy:
        X = X.copy() 
    m = np.median(X, axis=axis)
    q =  safe_asarray(mquantiles(X, prob=(0.25, 0.75), axis=axis))
    #iqr = (q[1, :] - q[0, :]) / 1.34898 # for Normal distributions, iqr = 1.34898*sigma
    iqr = np.sqrt((q[1, :] - q[0, :]))
    return (X - m) / iqr


def perfect_cluster():
    n_clusters = 10
    random_state = 42
    (Xn, X, zc, lc) = make_fabia_biclusters((100, 1000), n_clusters, (5, 25),
     (10, 210), 3, 0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
    model = FabiaBiclustering(n_clusters, 400, random_state=random_state)
    model.fit(X)
    Xr = np.dot(model.Z, model.L)
    f = plt.figure()
    a = plt.subplot(311)
    a.matshow(X, cmap=plt.cm.RdBu)
    a.set_title("noisless data")
    a = plt.subplot(312)
    a.matshow(Xn, cmap=plt.cm.RdBu)
    a.set_title("noisy data")
    a = plt.subplot(313)
    a.matshow(Xr, cmap=plt.cm.RdBu)
    a.set_title("reconstruction")
    plt.tight_layout(0)
    sklearn.metrics.consensus_score(model.biclusters_, (zc, lc))
    for a in f.axes: a.set_aspect('auto')
        
    Xn = np.loadtxt("/home/tom/Y.csv").T
    model = FabiaBiclustering(100, 500, random_state=random_state)
    model.fit(Xn)
   
    


def from_paper():
    n_clusters = 10
    random_state = 42
    n, m = 100, 1000
    (Xn, X, zc, lc) = make_fabia_biclusters((n, m), n_clusters, (5, 25),
     (10, 210), 3, 0.2, 2.0, 1.0, 0.2, 3.0, 1.0)
    model = FabiaBiclustering(n_clusters, 400, random_state=random_state, scale=False)
    model.fit(X)
    sklearn.metrics.consensus_score(model.biclusters_, (zc, lc))


    
    Xn = np.loadtxt("/home/tom/Y.csv").T
    model = FabiaBiclustering(100, 500, random_state=random_state)
    model.fit(Xn)
    
    L = np.loadtxt("/home/tom/L.csv").T
    



myfabia <- function(X,p=5,alpha=0.1,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0){
    l=ncol(X)
    n=nrow(X)
    rownames(X) <- rownames(X, do.NULL = FALSE, prefix = "gene")
    colnames(X) <- colnames(X, do.NULL = FALSE, prefix = "sample")
    rowna <- rownames(X)
    colna <- colnames(X)

    eps <- as.double(1e-3)
    eps1 <- as.double(1e-10)
    init_lapla <- 1.0
    init_psi <- 0.2
    iin <-  1.0/l

    XX <- as.vector(rep(1,n))
    if (center ==2) {
		cent <- apply(X, 1, median)
		X <- X - cent
	}
	else
		cent = 0
    X <- X - cent		
    XX <- as.vector(rep(1,n))
	if (norm == 2) {
		scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
		X <- scaleData*X	
	}
	else
		scaleData = 0
		
    #L <- matrix(random*rnorm(n*p),nrow=n,ncol=p)
	L <- as.matrix(t(read.table("/home/tom/L.csv")))

    lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)
    Psi <- init_psi*XX
	cyc <- as.integer(cyc)
	nL <- as.integer(nL)
	lL <- as.integer(lL)
	bL <- as.integer(bL)
	alpha <- as.double(alpha)
	non_negative <- as.integer(non_negative)
	p <- as.integer(p)
	spz <- as.double(spz)
	scale <- as.double(scale)
	lap <- as.double(lap)

	res <- .Call("fabic", X,Psi,L,lapla,cyc ,alpha,eps,eps1,spl,spz,scale,lap,nL,lL,bL,non_negative,PACKAGE="fabia")
 
    # INI call for biclusters
    vz <- iin*apply(res$E_SX_n,1,function(x) sum(x^2))
    vz <- sqrt(vz+1e-10)
    ivz <- 1/vz
    if(length(ivz)==1) {
        nZ <- ivz*res$E_SX_n
        noL <- vz*res$L
        res$lapla <- vz^2 * res$lapla
    }
    else {
        nZ <- ivz*res$E_SX_n
        noL <- t(vz*t(res$L))
        res$lapla <- sweep(res$lapla, 2, vz^2, "*")
    }
    return(new("Factorization", parameters=list("fabia",cyc,alpha,spl,spz,p,NULL,NULL,random,scale,norm,center,lap,nL,lL,bL,non_negative),
        n=n,p1=p,p2=p,l=l,center=cent,scaleData=scaleData,X=X,L=noL,Z=nZ,
        M=as.matrix(1),LZ=as.matrix(1),
        U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=res$Psi,lapla=res$lapla))
}
