
#include <hdf5.h>
#include <hdf5_hl.h>
#include <stdio.h>
#include <stdlib.h>

#include "util.h"


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
float* readMatrixFromHdf5(const char* fname, const char* dset, int* n, int* m) {

	hid_t fid = H5Fopen(fname, H5F_ACC_RDONLY, H5P_DEFAULT);
	
	if (fid < 0) {
		fprintf(stderr, "failed to open file\n");
		H5Eprint(H5E_DEFAULT, NULL);
		return 0;
	}
	
	int isPresent = H5LTfind_dataset(fid, dset);
	if (!isPresent) {
		H5Eprint(H5E_DEFAULT, NULL);
		return 0;
	}
		
	hsize_t dims[2];
	herr_t status = H5LTget_dataset_info(fid,dset,dims,NULL,NULL);
	if (status < 0) {
		H5Eprint(H5E_DEFAULT, NULL);
		return 0;
	}
 	*n = dims[0];
	*m = dims[1];
	
	float* data = (float*) malloc(dims[0] * dims[1] * sizeof(float));
	status = H5LTread_dataset_float(fid, dset, data);
	if (status < 0) {
		H5Eprint(H5E_DEFAULT, NULL);
		return 0;
	}
	
	H5Fclose(fid);
	return data;
}


/**
 * (Re)creates a HDF5-file, overwriting it if it already exists.
 */
void createHdf5File(const char* fname) {
	hid_t fid = H5Fcreate(fname, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	H5Fclose (fid);
}


/**
 * Stores a given matrix as dataset in a HDF5 file. Assumes the file was already created.
 * @param fname path to the HDF5 file
 * @param dset  name of the dataset
 * @param x     (flattened) matrix
 * @param n     size of the 1st dimension of x
 * @param m     size of the 2nd dimension of x
 * 
 * @return 0 on success, nonzero on error.
 */
int storeMatrixInHdf5(const char* fname, const char* dset, void* x, int n, int m) {

	hid_t fid = H5Fopen(fname, H5F_ACC_RDWR, H5P_DEFAULT);

	// we want to store the matrix in compressed format, so we can't just use H5LTmake_dataset
	//herr_t status = H5LTmake_dataset(fid, dset, int 2, dims, H5T_IEEE_F64LE, x);

	hsize_t dims[] = {n, m};
	hid_t dsid = H5Screate_simple (2, dims, NULL);
	hid_t plid  = H5Pcreate(H5P_DATASET_CREATE);
	hsize_t cdims[] = {n, m};
	int status = H5Pset_chunk (plid, 2, cdims);
	if (status < 0) {
		H5Eprint(H5E_DEFAULT, NULL);
		return 1;
	}
	
	status = H5Pset_deflate (plid, 7); // compression level 7
	if (status < 0) {
		H5Eprint(H5E_DEFAULT, NULL);
		return 1;
	}
	
	//hid_t memtype = isfloat ? H5T_IEEE_F64LE : H5T_IEEE_F32LE;
	hid_t memtype = H5T_IEEE_F32LE;
	hid_t ddid = H5Dcreate2(fid, dset, memtype, dsid, H5P_DEFAULT, plid, H5P_DEFAULT);
    status = H5Dwrite (ddid, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, x);
	if (status < 0) {
		H5Eprint(H5E_DEFAULT, NULL);
		return 1;
	}
	
    status = H5Sclose (dsid);
    status = H5Dclose (ddid);
    status = H5Pclose (plid);
    status = H5Fclose (fid);
	return 0;
}



void storeResults(const char* fname, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla) {
	
	/*
	 * delete the HDF5 file if it already exists.
	 * We need this because we can't delete datasets from already existing
	 * files. Thus, if we by accident try to store our results in an 
	 * already existing file, the program will fail to store its results
	 * there because a dataset of that name already exists.
	 * (there would be an option to edit existing datasets, but only if they're
	 * of the same size).
	 */
	createHdf5File(fname);
	storeMatrixInHdf5(fname, "Psi", Psi, n, 1);

	// NOTE: we get L, Z and lapla back in column major format
	float* tZ = toRowMajor(Z, k, l);
	storeMatrixInHdf5(fname, "Z", tZ, k, l);
	free(tZ);
	
	float* tL = toRowMajor(L, n, k);
	storeMatrixInHdf5(fname, "L", tL, n, k);
	free(tL);
	
	float* tlapla = toRowMajor(lapla, l, k);
	storeMatrixInHdf5(fname, "lapla", tlapla, l, k);
	free(tlapla);
}


void createRollbackPoint(int iter, float elapsedTime, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla)
{
	// create rollback point
	char fname[512];
	snprintf(fname, sizeof(fname), "fabia-rollback-k%d-cyc%d.hdf5", k, iter);
	storeResults(fname, k, n, l, L, Z, Psi, lapla);
}


/*
int main(int argc, const char** argv) {
	int n, m;
	float* x = readMatrixFromHdf5("/media/scratch/data/cifar10bw.hdf5", "testx", &n, &m);
	printf("first: %f %f %f\n", x[0], x[1], x[2]);
	printf("dim: %d %d\n", n, m);
	
	int status = storeMatrixInHdf5("test.hdf5", "testx", x, n, m);
	printf("status: %d\n", status);
	
	return 0;
}
*/
