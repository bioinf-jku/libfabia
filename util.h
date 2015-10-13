#ifndef UTIL_H
#define UTIL_H

struct timespec;

float calculateElapsedTime(struct timespec t0);

 // uniform random in (0, 1]
double rand_unif(void);

// random from a N(0, 1) Gaussian
double rand_normal(void);

/// Transforms a matrix from Row-Major to Column-Major.
float* toColumnMajor(float* a, int n, int m);

/// Transforms a matrix from Column-Major to Row-Major.
float* toRowMajor(float* a, int n, int m);


/// Copies a matrix from Column-Major to Row-Major.
void copyRowMajor(float* src, float* dest, int n, int m);

/// Transposes a matrix in Column-Major format.
float* transposeMatrixCM(float* a, int n, int m);

/// Prints a column major matrix.
void printMatrixCM(float* a, int n, int m);

/// Prints a row-major matrix
void printMatrixRM(float* a, int n, int m);

float* asFloat(double* a, int n);


double* asDouble(float* a, int n);

#endif
