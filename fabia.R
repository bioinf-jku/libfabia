runFabia <- function(X, p=5, cyc=500, alpha=0.01, spl=0, spz=0.5, scale=1.0, nThreads=1) {
    dyn.load("libfabia.so")
	l=ncol(X)
	n=nrow(X)
	eps <-1e-3
	eps1 <- 1e-10
    lap = 1.0
	init_lapla <- 1.0
	init_psi <- 0.2
	
	XX <- as.vector(rep(1,n))
	L <- matrix(rnorm(n*p),nrow=n,ncol=p)
    Z <- matrix(rnorm(l*p),nrow=p,ncol=l)
	lapla <-  init_lapla*matrix(1,nrow=l,ncol=p)
	Psi <- init_psi*XX
	Psi <- init_psi*XX
	
    system.time(resC <- .C("fabia_cm_d_r", as.integer(p), as.integer(n), as.integer(l),
    X=as.double(X), Psi=as.double(Psi), L=as.double(L), Z=as.double(Z),
    lapla=as.double(lapla), cyc, as.double(alpha), as.double(eps), as.double(spl),
    as.double(spz), as.integer(scale), as.double(lap), as.integer(10), nThreads))

    # resizing results (dimensions get lost during .C)
    res=list()
    res$L <- matrix(resC$L, nrow=n, ncol=p)
    res$Z <- matrix(resC$Z, nrow=p, ncol=l)
    res$lapla <- matrix(resC$lapla, nrow=l,ncol=p)
    res$Psi <- resC$Psi
	res$p <- p

	# rescaling Z so it has variance 1
	vz <- apply(res$Z,1,function(x) sum(x^2)) / ncol(X)
	vz <- sqrt(vz+1e-10)
	ivz <- 1.0/vz
	if(length(ivz)==1) {
		res$Z <- ivz*res$Z
		res$L <- vz*res$L
		
	}else {
		res$Z <- ivz*res$Z
		res$L <- sweep(res$L, 2, vz, "*")
		res$lapla <- sweep(res$lapla, 2, vz^2, "*")
	}
    dyn.unload("libfabia.so")
	return(res)
}

#runFabia( matrix(rnorm(20), ncol=4))
