#include <stdlib.h>
#include <stdio.h>  
#include <limits.h> 
#include <math.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <omp.h>


#include "util.h"
#include "hdf5tools.h"


#include "f2c.h"
#define integer int
#include "clapack.h"
#include "blas.h"

// change the "s" to "d" when switching from float to double
#define F77(x) s ## x ## _


static const float MACHINE_EPS = 1e-7;

static inline void invertCholesky(float* a, int n) { 
	int info;
	F77(potrf)("l", &n, a, &n, &info);
	assert(!info);
	F77(potri)("l", &n, a, &n, &info);
	assert(!info);
	// dpotri only solves the lower triangle of the symmetric matrix
	for (int i1 = 0; i1 < n; ++i1) {
		for (int i2 = i1+1; i2 < n; ++i2) {
			//a[i2 + i1*n] = a[i1 + i2*n];
			a[i1 + i2*n] = a[i2 + i1*n];
		}
	}
}


extern void updateUI(int iter, float elapsedTime, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla);
extern float calculateElapsedTime(struct timespec t0);


// needed for BLAS/LAPACK calls
static const float ONE = 1;
static const float ZERO = 0;
static const int ONE_INT = 1;
	

/**
 * Updates the estimate E(z|x) for a given datapoint x, using the current LPsi/LPsiL matrix.
 * Uses the work-matrices tLPsiL and iLPsiLX, whose contents will be destroyed by this call.
 * 
 * If sum2 != NULL, then the estimates for lapla, as well as sum2 will also be 
 * calculated.
 */
static inline void estimateZ(int n, int k, const float* x, float* z, float* lapla, float* LPsi, float* LPsiL, float* tLPsiL, float* iLPsiLX,
	float* sum2, float spz, float lap)
{
	// x = LPsiL+diag(lapla[j,])
	for (int i = 0; i < k; ++i) {
		memcpy(tLPsiL+i*k, LPsiL+i*k, k*sizeof(float));
		tLPsiL[i*k + i]+= lapla[i];
	}

	invertCholesky(tLPsiL, k); // tmp <- chol2inv(chol(x))
	F77(symm)("l", "l", &k,&n,&ONE,tLPsiL,&k, LPsi,&k,&ZERO,iLPsiLX,&k); //z <- (tmp %*% t(LPsi))%*% x
	F77(gemv)("n", &k, &n, &ONE, iLPsiLX, &k, x, &ONE_INT, &ZERO,  z, &ONE_INT);
	
	if (!sum2)
		return;
			
	// note: from here on out, we use tLPsiL to store zz, so we don't have to copy it
	F77(ger)(&k, &k, &ONE, z, &ONE_INT, z, &ONE_INT, tLPsiL, &k); //zz <-  tmp + z %*% t(z)
	for (int i1 = 0; i1 < k; ++i1) { // sum2 <- sum2 + zz
		for (int i2 = 0; i2 < k; ++i2) {
			sum2[i1*k +i2] += tLPsiL[i1*k + i2];
		}
	}
			
	for (int i1=0;i1<k;i1++) { //laj <- (epsv+diag(zz))^(-spz); laj[which(laj<lap)] <- lap; lapla[j,] <- lajr
		float s = pow(MACHINE_EPS+tLPsiL[i1*k+i1],-spz);
		if (s<lap) {
			lapla[i1] = lap;
		} else {
			lapla[i1] = s; 
		}
	}
}


/**
 * Runs the main FABIA algorithm. Expects all  matrices in column-major layout
 * and operates on float matrices.
 * 
 * @param X	matrix of n * l, with nn datapoints in its columns
 * @param Psi	vector of n 
 * @param L	matrix of n * k
 * @param Z	matrix of k * l
 * @param lapla	matrix of k * l   (NOTE: this is the transpose of what it is in R!)
 * @param cyc	the number of EM cycles
 * @param alpha	parameter for the Laplace prior
 * @param eps	epsilon used for regularization
 * @param spl	parameter to tune extra-sparseness of l
 * @param spz	parameter to tune extra-sparseness of z
 * @param scale	scale parameter, as in original FABIA
 * @param lap	minimal value of the variational parameter
 * @param verbose  if non-zero, print status information every $verbose iterations.
 * @param nthreads the number of threads that will be used
 */
void fabia_cm_f(const int k, const int n, const int l, 
	const float* restrict X, float* restrict Psi, float* restrict L, float* restrict Z, float* restrict lapla,
	int cyc, float alpha, float eps, float spl, float spz, int scale, float lap, int verbose, int nthreads)
{
	float *XZ = (float*) malloc(n * k * sizeof(float));
	float *sum2_ = (float*) malloc(nthreads * k * k * sizeof(float));
	float *XX = (float*) malloc(n * sizeof(float)); 
	float *tLPsiL_ = (float*) malloc(nthreads * k * k * sizeof(float));
	float *LPsiL = (float*) malloc(k * k * sizeof(float));
	float *LPsi = (float*) malloc(k * n * sizeof(float)); 
	float *iLPsiLX_ = (float*) malloc(nthreads * n * k * sizeof(float));
    
	// Note: LPsi is transposed wrt Sepp's implementation!
	// This is because then we can make use of the fact that inv(LPsiL) is
	// symmetric when multiplying it with LPsi in the inner (j) loop	
	
	if (!iLPsiLX_) { 	// if the last allocation failed, we ran out of memory
		fprintf(stderr, "Out of memory\n");
		free(tLPsiL_);
		free(LPsi);
		free(LPsiL);
		free(XX);
		free(sum2_);
		free(XZ);
		memset(L, 0, n*k*sizeof(float));
		memset(Z, 0, k*l*sizeof(float));
		return;
	}

	// these were parameters in the original FABIA, but we never use them.
	//const int nL = 0;
	//const int lL = 0;
	//const int bL = 0;
	const int non_negative = 0;
	
	omp_set_num_threads(nthreads);
	
	struct timespec t0;
	clock_gettime(CLOCK_MONOTONIC, &t0);


	if (lap<eps)
		lap = eps;
	
	memset(XX, 0, n*sizeof(float));
	for (int i1 = 0; i1 < l; ++i1) {
		for (int i2 = 0; i2 < n; ++i2)
			XX[i2] += X[i2 + i1*n] * X[i2 + i1*n];
	}
	for (int i2 = 0; i2 < n; ++i2) {
		XX[i2] /= l;
        XX[i2] = max(XX[i2], eps);
    }
    
	float t = 0.0;
	for (int iter = 1; iter <= cyc; ++iter) {

		for (int i1 = 0; i1 < k; i1++) { //LPsi<-diag(1/Psi)%*%L
			for (int i2=0;i2<n;i2++) {
				LPsi[i2*k + i1] =  L[i2 + i1*n]/Psi[i2]; // LPsi gets transposed here
			}
		}

 		F77(gemm)("t","t",&k,&k,&n, &ONE,L, &n, LPsi,&k,&ZERO,LPsiL,&k); //LPsiL = t(L) %*% LPsi	
 		
 		// zero out all sums (TODO: maybe faster if ran in parallel?)
		memset(sum2_, 0, nthreads*k*k*sizeof(float));
		
		// the first k*k elements of sum2 will be used to aggregate, so initialize diagonal
		for (int i=0; i < k; ++i)
			sum2_[i*k + i] = eps;
		
		#pragma omp parallel for
		for (int j = 0; j < l; ++j) {
			int tid = omp_get_thread_num();
			float* tLPsiL = tLPsiL_ + (tid*k*k);
			float* iLPsiLX = iLPsiLX_ + (tid*n*k);
			float* sum2 = sum2_ + (tid*k*k);
			estimateZ(n, k, X+n*j, Z+j*k, lapla+j*k, LPsi, LPsiL, tLPsiL,iLPsiLX, sum2, spz, lap);
		}	
		// aggregate/reduce sum2 into their first sub-arrays
		#pragma omp parallel for
		for(int i = 0; i < k*k; ++i)
			for (int tid = 1; tid < nthreads; ++tid)
				sum2_[i] += sum2_[i + tid*k*k];
						
		invertCholesky(sum2_, k);  // sll <- chol2inv(chol(sum2))
		F77(gemm)("n","t",&n,&k,&l,&ONE,X,&n, Z,&k,&ZERO,XZ,&n);
		F77(gemm)("n","n",&n,&k,&k,&ONE,XZ,&n, sum2_,&k,&ZERO,L,&n); // L <- sum1%*%sll
		
		//L <- L - alpha*Psi*sign(L)*abs(nk_one*eps+L)^{-spl}
		for (int i2 = 0; i2 < k; ++i2) {
			for (int i1 = 0; i1 < n; ++i1) {
				float s = L[i1+ n*i2];
				float sgn = (s > 0) - (s < 0);
				if ((non_negative<=0) || (sgn>0)) {
					float t = fabs(Psi[i1]*alpha*pow(MACHINE_EPS+fabs(s), -spl));
					if (fabs(s) > t)
						L[i1+ n*i2] -= sgn * t;
					else
						L[i1+ n*i2] = 0.0;
				} else
					L[i1+ n*i2] = 0.0;
			}
		}
		
		// TODO: deal with nL here! (check in the original code ~ line 470)
		// TODO: deal with lL here! (check in the original code ~ line 504)
		
		//Psi <- epsn+abs(XX - diag(tcrossprod(L,sum1)/n))
		t = 0.0;
		for (int i1 = 0; i1 < n ;++i1) {
			float s = 0.0;
			for (int i2 = 0; i2 < k; ++i2) {
				s += L[i1 + i2*n] * XZ[i1 + i2*n];
			}
			if (fabs(s)>t)
				t = fabs(s);
			Psi[i1] = XX[i1] - s/l;
			if (Psi[i1] < eps) {
				Psi[i1] = eps;
			}
		}
		if (t < eps ) {
				
			for (int i1=0;i1<n;i1++) {
				Psi[i1] = eps;
			}
			for (int j=0; j < l; j++) {
				for (int i1=0; i1 < k; i1++)
					lapla[j*k+i1] = eps; 
			}
			printf("Last update was %f, which is smaller than %f, so I'm bailing out\n", t, eps);
			break;
		}
		
		if (scale) {
			for (int i = 0; i < k; i++){
				float s = 0.0;
				for (int j=0;j<n;j++) {
					s+=L[i*n + j]*L[i*n + j];
				}
				s= 1.0 / (sqrt(s/n)+MACHINE_EPS);
				for (int j=0;j<n;j++) {
					L[i*n + j]*=s;
				}
				s = pow(s*s,-spz);
				for (int j=0;j<l;j++) {
					lapla[j*k + i]*=s;
				}
			}
		}
		
		// check (and re-initialize) "all-zeroes" bicluster
		for (int i = 0; i < k; ++i) {
			int isZero = 1;
			for (int j = 0; j < n; ++j) {
				if (L[j + i*n] != 0) {
					isZero = 0;
					break;
				}
			}
			
			if (isZero) {
				for (int j = 0; j < n; ++j)
					L[i*n + j] = (float) rand_normal();
				for (int j = 0; j < l; ++j)
					lapla[j*k + i] = 1.0;
			}
		}		
		
		if (verbose && (iter % verbose == 0)) {
			updateUI(iter, calculateElapsedTime(t0), k, n, l, L, Z, Psi, lapla);
		}
	}
	
	// last update
	if (t >= eps ) {
		for (int i1 = 0; i1 < k; i1++) { 		//LPsi<-diag(1/Psi)%*%L
			for (int i2=0;i2<n;i2++) {
				LPsi[i2*k + i1] =  L[i2 + i1*n]/Psi[i2]; // LPsi gets transposed here
			}
		}
		F77(gemm)("t","t",&k,&k,&n, &ONE,L, &n, LPsi,&k,&ZERO,LPsiL,&k); //LPsiL = t(L) %*% LPsi	

		#pragma omp parallel for
		for (int j = 0; j < l; ++j) {
			int tid = omp_get_thread_num();
			float* tLPsiL = tLPsiL_ + (tid*k*k);
			float* iLPsiLX = iLPsiLX_ + (tid*n*k);

			estimateZ(n, k, X+n*j, Z+j*k, lapla+j*k, LPsi, LPsiL, tLPsiL,iLPsiLX, NULL, spz, lap);
		}
	}
	else
		memset(Z, 0, k*l*sizeof(float));


	free(XX); 
	free(XZ);
	free(sum2_);
	free(iLPsiLX_);
	free(tLPsiL_);
	free(LPsiL);
	free(LPsi);
	return;
}
