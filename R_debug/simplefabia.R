myfabi <- function(X,p=5,alpha=0.01,cyc=500,spl=0,spz=0.5,non_negative=0,random=1.0,center=2,norm=1,scale=0.0,lap=1.0,nL=0,lL=0,bL=0){
        l=ncol(X)
        n=nrow(X)
        
        if (center ==2) {
			cent <- apply(X, 1, median)
			X <- X - cent
		}
		
        XX <- as.vector(rep(1,n))
		if (norm == 2) {
			scaleData <-  1/sqrt(apply(X,1,function(x) {quantile(x,0.75) - quantile(x,0.25)})+0.001*XX)
			X <- scaleData*X	
		}


        eps <- as.double(1e-3)
        #eps1 <- as.double(1e-10)
        eps1 <- 2.2204460492503131e-15
        init_lapla <- 1.0
        init_psi <- 0.2
        iin <-  1.0/l
        lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)
        Psi <- init_psi*XX

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
		XX <- apply(X*X, 1, sum) / as.numeric(l)
		#XX <- apply(X, 1, var)
		#L <- matrix(rnorm(n*p),nrow=n,ncol=p)
		L <- as.matrix(t(read.table("/home/tom/L.csv")))
        #L <- (0.5*XX^(.5))*nk_one

		for (i in 1:cyc){
			LPsi<-diag(1/Psi)%*%L
			LPsiL<-crossprod(L,LPsi)
			sum1<- nk_zero
			sum2<- eps*kk_one
			for (j in 1:l){
				#tmp <- chol2inv(chol(LPsiL+diag(lapla[j,])))
				tmp <- solve(LPsiL+diag(lapla[j,]))
				x_j <- as.vector(X[,j])
				E_SX_n[,j] <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
				e_ssxx_n <-  tmp + tcrossprod(E_SX_n[,j],E_SX_n[,j])
				sum1 <- sum1 +  tcrossprod(x_j,E_SX_n[,j])
				sum2 <- sum2 + e_ssxx_n
				laj <- (epsv+diag(e_ssxx_n))^(-spz)
				laj[which(laj<lap)] <- lap
				lapla[j,] <- laj
			}
			#sll <- chol2inv(chol(sum2))
			sll <- solve(sum2)
			L <- sum1%*%sll
			ddL <- alpha*Psi*sign(L)*abs(nk_one*eps+L)^{-spl}
			#ddL <- alpha*Psi*L
			L <- L - ddL
			L[which(abs(L)<abs(ddL))] <- 0
			Psi <- abs(XX - diag(tcrossprod(L,sum1)/l))
			Psi[Psi < eps] = eps
			if (i %% 20==0) { print(i)}
			if (FALSE) {
			print(i)
			print(cat("XX sum:  ", sum(abs(XX)), " "))
			print(cat("xz sum:  ", sum(abs(sum1)), " "))
			print(cat("sll sum: ", sum(abs(sll)), " "))
			print(cat("sum2 sum:", sum(abs(sum2)), " "))
			print(cat("z sum:   ", sum(abs(E_SX_n)), " "))
			print(cat("l sum:   ", sum(abs(L)), " "))
			print(cat("psi sum: ", sum(abs(Psi)), " "))
			print(cat("lap sum: ", sum(abs(lapla)), " "))
			}
	}

#print("final L:")
#print(L[1:5, 1:5])
#print("final lapla:")
#print(lapla[1:5, 1:5])
#print("final Z:")
#print(E_SX_n[1:5, 1:5])

    LPsi<-diag(1/Psi)%*%L
    LPsiL<-crossprod(L,LPsi)
    for (j in 1:l){
		#tmp <- chol2inv(chol(LPsiL+diag(lapla[j,])))
		tmp <- solve(LPsiL+diag(lapla[j,]))
        x_j <- as.vector(X[,j])
        e_sx_n <- as.vector(tcrossprod(tmp,LPsi)%*%x_j)
        E_SX_n[,j] <- e_sx_n
        E_SSXX_n[[j]] <-  tmp + tcrossprod(e_sx_n,e_sx_n)
    }

    # INI call for biclusters
    vz <- apply(E_SX_n,1,function(x) sum(x^2)) / l
    vz <- sqrt(vz+1e-10)
    Z <- E_SX_n
    if(length(vz)==1) {
        Z <- E_SX_n / vz
        L <- vz*L
        lapla <- vz^2 * lapla
     }
    else {
        Z <- E_SX_n / vz
        L <- t(vz*t(L))
        lapla <- sweep(lapla, 2, vz^2, "*")
     }
	#return (list(L=L, Z=Z, Psi=Psi, lapla=lapla))
    return(new("Factorization", parameters=list("fabia",cyc,alpha,spl,spz,p,NULL,NULL,random,scale,norm,center,lap,nL,lL,bL,non_negative),
        n=n,p1=p,p2=p,l=l,center=0,scaleData=0,X=X,L=L,Z=Z,
        M=as.matrix(1),LZ=as.matrix(1),
        U=as.matrix(1),avini=as.vector(1),xavini=as.vector(1),ini=as.matrix(1),Psi=Psi,lapla=lapla))
}


#install.packages("clue")


consensus_score <- function(n, l, colsA, rowsA, colsB, rowsB) {

    convert <- function(n, l, sampleList, featureList) {
        p = length(sampleList)
        m_rows = matrix(0,n,p)
        m_cols = matrix(0,l,p)
        for (i in seq_along(sampleList)) {
            m_cols[sampleList[[i]], i] <- 1
            m_rows[featureList[[i]], i] <- 1
        }
        return (list(m_rows, m_cols))
    }
    
    tmpA = convert(n, l, colsA, rowsA)
    rowsA = tmpA[[1]]
    colsA = tmpA[[2]]
    
    tmpB = convert(n, l, colsB, rowsB)
    rowsB = tmpB[[1]]
    colsB = tmpB[[2]]

    jaccard <- function(a_rows, a_cols, b_rows, b_cols) {
        is = sum(a_rows * b_rows) * sum(a_cols * b_cols)
        as = sum(a_rows) * sum(a_cols)
        bs = sum(b_rows) * sum(b_cols)
        return (is / (as+bs - is))
    }
    
    # jacarrad consenus matrix:
    pA = ncol(colsA)
    pB = ncol(colsB)
    jamat <- matrix(0,pA,pB)
    for (i in 1:pA) {
        for (j in 1:pB) {
            jamat[i, j] = jaccard(rowsA[, i], colsA[, i], rowsB[, j], colsB[, j])
        }
    }
   
    require(clue)
    # solve_LSAP needs #colums > #rows
    if (ncol(jamat) < nrow(jamat))
        jamat = t(jamat)
    y = solve_LSAP(jamat, maximum=TRUE)
    score = sum(jamat[cbind(seq_along(y), y)])
    return (list(score / max(pA, pB), jamat))
}






library(fabia, verbose=F, quietly=T)
dat <- makeFabiaData(n = 1000,l= 100,p = 10,f1 = 5,f2 = 5,of1 = 5,of2 = 10,  sd_noise = 3.0,sd_z_noise = 0.2,mean_z = 2.0,sd_z = 1.0, sd_l_noise = 0.2,mean_l = 3.0,sd_l = 1.0)

load("/home/tom/fabdata.rdata")
X = dat$X
#X = (X - rowMeans(X)) 
#X = X / apply(X, 1, sd)
res <- myfabi(X,p=13, alpha=0.01, cyc=500, norm=0, center=0)

# true results
zz = matrix(0, nrow=10, 100)
for (i in seq_along(dat$ZC))
    zz[i, dat$ZC[[i]]] = 1
ll = matrix(0, nrow=10, 1000)
for (i in seq_along(dat$LC))
    ll[i, dat$LC[[i]]] = 1   

write.table(dat$X, file="X.csv", , row.names=F, col.names=F)
write.table(zz, file="/home/tom/row_orig.csv", row.names=F, col.names=F)
write.table(ll, file="/home/tom/col_orig.csv", row.names=F, col.names=F)

# fabia results
write.table(res@L, file="/home/tom/Lres.csv", row.names=F, col.names=F)
write.table(res@Z, file="/home/tom/Zres.csv", row.names=F, col.names=F)


library(fabia)
bic = extractBic(res)
sampleList = list()
featureList = list()
for (i in 1:bic$np) {
    sampleList[[i]] = bic$numn[i,]$numnp
    featureList[[i]] = bic$numn[i,]$numng
}

n=1000
l=100
s = consensus_score(n, l, sampleList, featureList, dat$ZC, dat$LC)
print(s)

