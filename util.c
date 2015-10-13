#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "hdf5tools.h"
#include "util.h"

#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

double rand_unif(void) {return (rand()+1.0)/(RAND_MAX+1.0);} // random in (0, 1]

// generates random samples from a 0/1 Gaussian via Box-Mueller
double rand_normal(void) {return sqrt(-2.0*log(rand_unif())) * cos(2.0*M_PI*rand_unif());} 


/// Transforms a matrix from Row-Major to Column-Major.
float* toColumnMajor(float* a, int n, int m) {
	float* b = (float*) malloc(n*m*sizeof(float));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j)
			b[j*n + i] = a[i*m +j];
	}
	return b;
}


/// Copies a matrix from Column-Major to Row-Major.
void copyRowMajor(float* src, float* dest, int n, int m) {
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j)
			dest[i*m +j] = src[j*n + i];
	}
}


/// Transforms a matrix from Column-Major to Row-Major.
float* toRowMajor(float* a, int n, int m) {
	float* b = (float*) malloc(n*m*sizeof(float));
	copyRowMajor(a, b, n, m);
	return b;
}


/// Transposes a matrix in Column-Major format.
float* transposeMatrixCM(float* a, int n, int m) {
	float* b = (float*) malloc(n*m*sizeof(float));
	for (int i = 0; i < n; ++i) {
		for (int j = 0; j < m; ++j)
			b[m*i + j] = a[j*n +i];
	}
	return b;
}


/// Prints a column major matrix.
void printMatrixCM(float* a, int n, int m) {
	for (int i = 0; i < n; ++i) {
		for (int j =0 ; j < m; ++j) 
			printf("%1.3f ", a[i + j*n]);
		printf("\n");
	}
}


/// Prints a row-major matrix
void printMatrixRM(float* a, int n, int m) {
	for (int i = 0; i < n; ++i) {
		for (int j =0 ; j < m; ++j) 
			printf("%1.3f ", a[i*m + j]);
		printf("\n");
	}
}


float* asFloat(double* a, int n) {
	float* retval = (float*)malloc(n * sizeof(float));
	for (size_t i = 0; i < n; ++i)
		retval[i] = (float) a[i];
	return retval;
}


double* asDouble(float* a, int n) {
	double* retval = (double*)malloc(n * sizeof(double));
	for (size_t i = 0; i < n; ++i)
		retval[i] = (double) a[i];
	return retval;
}


float calculateElapsedTime(struct timespec t0) {
	struct timespec t;
	clock_gettime(CLOCK_MONOTONIC, &t);
	float elapsed = (t.tv_sec - t0.tv_sec);
	elapsed += (t.tv_nsec - t0.tv_nsec) / 1000000000.0;	
	return elapsed;
}
