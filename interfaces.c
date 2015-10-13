

#include "fabia.h"
#include "util.h"
#include <stdlib.h>


/// This is the function callable through the .C interface in R
/// NOTE: * internally, the function will still run in single-precision
///       * lapla is transposed wrt. the original implementation!
void fabia_cm_d_r(int* k, int* n, int* l, double* X, double* Psi, double* L, double* Z, double* lapla,
	int* cyc, double* alpha, double* eps, double* spl, double* spz, int* scale, double* lap, int* verbose, int* nthreads)
{
	float* Xf = asFloat(X, n[0]*l[0]);
	float* Psif = asFloat(Psi, n[0]);
	float* Lf = asFloat(L, n[0]*k[0]);
	float* Zf = asFloat(Z, k[0]*l[0]);
	float* laplaf = asFloat(lapla, k[0]*l[0]);
	
	fabia_cm_f(k[0], n[0], l[0], Xf, Psif, Lf, Zf, laplaf, cyc[0], (float) alpha[0],
		(float)eps[0], (float)spl[0], (float)spz[0], scale[0], (float) lap[0], verbose[0], nthreads[0]);
	
	for (int i = 0; i < n[0]*l[0]; ++i) X[i] = (double) Xf[i];
	for (int i = 0; i < n[0]; ++i) Psi[i] = (double) Psif[i];
	for (int i = 0; i < n[0]*k[0]; ++i) L[i] = (double) Lf[i];
	for (int i = 0; i < k[0]*l[0]; ++i) Z[i] = (double) Zf[i];
	for (int i = 0; i < k[0]*l[0]; ++i) lapla[i] = (double) laplaf[i];
		
	free(Xf);
	free(Psif);
	free(Lf);
	free(Zf);
	free(laplaf); 
}


/// interface that accepts row-major data
void fabia_rm_f(const int k, const int n, const int l, 
	float* restrict X, float* restrict Psi, float* restrict L, float* restrict Z, float* restrict lapla,
	int cyc, float alpha, float eps, float spl, float spz, int scale, float lap, int verbose, int nthreads)
{
	float* Xf = toColumnMajor(X, n, l);
	float* Lf = toColumnMajor(L, n, k);
	float* Zf = toColumnMajor(Z, k, l);
	float* laplaf = toColumnMajor(lapla, k, l);
	
	fabia_cm_f(k, n , l, Xf, Psi, Lf, Zf, laplaf, cyc, (float) alpha,
		(float)eps, (float)spl, (float) spz, scale, (float) lap, verbose, nthreads);
		
	copyRowMajor(Xf, X, n, l);
	copyRowMajor(Lf, L, n, k);
	copyRowMajor(Zf, Z, k, l);
	copyRowMajor(laplaf, lapla, k, l);
		
	free(Xf);
	free(Lf);
	free(Zf);
	free(laplaf);
}


