#ifndef FABIA_H
#define FABIA_H


/**
 * Runs the main FABIA algorithm. Expects all  matrices in column-major layout
 * and operates on float matrices.
 * 
 * @param X	matrix of n * l, with l datapoints in its columns
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
 * @param verbose  if non-zero, print status information
 * @param nthreads the number of threads that will be used
 */
void fabia_cm_f(const int k, const int n, const int l, 
	const float* X, float*  Psi, float*  L, float*  Z, float*  lapla,
	int cyc, float alpha, float eps, float spl, float spz, int scale, float lap, int verbose, int nthreads);
	

void fabia_cm_d(const int k, const int n, const int l, 
	double* restrict X, double* restrict Psi, double* restrict L, double* restrict Z, double* restrict lapla,
	int cyc, double alpha, double eps, double spl, double spz, int scale, double lap, int verbose, int nthreads);
	
/// Approximate 
void approx_fabia_cm_f(const int k, const int n, const int l, 
	float*  X, float*  Psi, float*  L, float*  Z, float*  lapla,
	int cyc, float alpha, float eps, float spl, float spz, int scale, float lap, int verbose, int nthreads);
#endif
