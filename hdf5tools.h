#ifndef MY_HDF5TOOLS_H
#define MY_HDF5TOOLS_H

#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>


/**
 * Reads a matrix from a HDF5 file.
 * It is assumed that the matrix is actually contained in the file, that
 * it's datatype is float and that it is 2dimensional.
 * Note that the data will be converted to float when returning it!
 * @param fname path to the HDF5 file
 * @param deset name of the dataset in the file.
 * @param n     returns the size of the 1st dimension of the matrix
 * @param m     returns the size of the 2nd dimension of the matrix
 *
 * @return flattened array containing the data. The user is responsible for
 *         free()-ing the data.
 */
float* readMatrixFromHdf5(const char* fname, const char* dset, int* n, int* m);


/**
 * (Re)creates a HDF5-file, overwriting it if it already exists.
 */
int createHdf5File(const char* fname) ;


/**
 * Stores a given matrix as dataset in a HDF5 file.
 * @param fname path to the HDF5 file
 * @param dset  name of the dataset
 * @param x     (flattened) matrix
 * @param n     size of the 1st dimension of x
 * @param m     size of the 2nd dimension of x
 *
 * @return 0 on success, nonzero on error.
 */
int storeMatrixInHdf5(const char* fname, const char* dset, void* x, int n, int m);

void storeResults(const char* fname, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla);

void createRollbackPoint(int iter, float elapsedTime, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla);

#endif
