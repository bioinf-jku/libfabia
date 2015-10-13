#include "fabia.h"
#include "hdf5tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <getopt.h>
#include <time.h>
#include "util.h"


void updateUI(int iter, float elapsedTime, int k, int n, int l, float* L, float* Z, float* Psi, float* lapla) {
	printf("Elapsed time: %7.1fs Cycle: %d\n", elapsedTime, iter);
	createRollbackPoint(iter, elapsedTime, k, n, l, L, Z, Psi, lapla);
}


/// samples a given fraction of datapoints from matrix a, destroying a in the process.
float* sample(float* a, int* n, int m, float fraction) {
	if (fraction < 1.0) {
		int* idx = (int*)malloc(*n*sizeof(int));
		int cnt = 0;
		for (int i = 0; i < *n; ++i) {
			idx[i] =  rand_unif() < fraction;
			cnt += idx[i];
		}
		
		float* b = (float*) malloc(cnt * m * sizeof(float));
		int writeidx = 0;
		for (int i = 0; i < *n; ++i) {
			if (!idx[i])
				continue;
			memcpy(b+writeidx*m, a + i*m, m*sizeof(float));
			++writeidx;
		}
		free(idx);		
		free(a);
		*n = cnt;
		return b;
	} else
		return a;
}


float* readMatrixFromRestartFile(const char* fname, const char* dset, int n, int m) {
	int sizeA, sizeB;
	float* tA = readMatrixFromHdf5(fname, dset, &sizeA, &sizeB);
	if (!tA || sizeA != n || sizeB != m) {
		fprintf(stderr, "Invalid restart file, unable to read %s\n", dset);
		return 0;
	}
	float *tmp = toColumnMajor(tA, n, m);
	free(tA);
	return tmp;
}


const char HELP_MESSAGE[] = 
"Runs FABIA on a dataset stored in HDF5, stores result in HDF5 file.\n\n"
" USAGE: fabia inputfile dsetname outputfile\n\n"
" parameters: \n"
" -v/--verbose n     print progress status and create rollback point every n cycles .\n"
" -h/--help 	     show help message.\n"
" -k/--nbicluster k  look for n biclusters.\n"
" -n/--nthreads nt   use a total of nt threads.\n"
" -c/--cyc  c        run FABIA for c cycles.\n"
" -a/--alpha  a      set Laplace prior to a.\n"
" -s/--sample f      sample fraction f from the data, instead of using all of it.\n"
" -r/--restart fn    restart computation using results from file fn.\n"
" -t/--seed s        Seed for the random number generator.\n"
" -l/--scaleL         If set, the L matrix will be rescaled after each iteration.\n"
" -x/--approx         If set, use an approximate algorithm.\n"
"\n";
	 
int main(int argc, char* const* argv) {
	
	int verbose = 0;
	const char* inputfile = NULL;
	const char* outputfile = NULL;
	const char* dataset = NULL;
	const char* restartfile = NULL;
	int k = 10;
	int nthreads = 1;
	int cyc = 1000;
	float alpha = 0.01;
	float sampleFraction = 1.0;
	unsigned int seed = time(NULL);
	int scale = 0; // wether to scale L values (this should prevent biclusters from dying out
	int approximate = 0;
	
	while (1) {
		static struct option long_options[] = {
			{"verbose",    required_argument, 0, 'v'},
			{"help",       no_argument,       0, 'h'},
			{"nbicluster", required_argument, 0, 'k'},
			{"nthreads",   required_argument, 0, 'n'},
			{"cyc",        required_argument, 0, 'c'},
			{"alpha",      required_argument, 0, 'a'},
			{"sample",     required_argument, 0, 's'},
			{"restart",    required_argument, 0, 'r'},
			{"seed",       required_argument, 0, 't'},
			{"scaleL",     no_argument,       0, 'l'},
			{"approx",     no_argument,       0, 'x'},
			{0, 0, 0, 0}
		};

		int opt_idx = 0;
		int c = getopt_long (argc, argv, "v:f:d:k:n:c:a:ho:s:r:t:lx", long_options, &opt_idx);

		if (c == -1)
			break;
		else if (c == 'v' || (opt_idx && !strcmp("verbose", long_options[opt_idx].name)))
			verbose = atoi(optarg);
		else if (c == 'k' || (opt_idx && !strcmp("nbicluster", long_options[opt_idx].name))) {
			k = atoi(optarg);
			if (k < 2) {
				fprintf(stderr, "Invalid value for k\n");
				return -1;
			}
		} else if (c == 'n' || (opt_idx && !strcmp("nthreads", long_options[opt_idx].name))) {
			nthreads = atoi(optarg);
			if (nthreads < 1)  {
				fprintf(stderr, "Invalid value for nthreads\n");
				return -1;
			}
		} else if (c == 'c' || (opt_idx && !strcmp("cyc", long_options[opt_idx].name))) {
			cyc = atoi(optarg);
			if (cyc < 1)  {
				fprintf(stderr, "Invalid value for cyc\n");
				return -1;
			}
		} else if (c == 'a' || (opt_idx && !strcmp("alpha", long_options[opt_idx].name))) {
			alpha = atof(optarg);
			if (alpha < 1e-8) {
				fprintf(stderr, "Invalid value for alpha\n");
				return -1;
			}
		} else if (c == 's' || (opt_idx && !strcmp("sample", long_options[opt_idx].name))) {
			sampleFraction = atof(optarg);
			if (1e-8 > sampleFraction || sampleFraction > 1.0) {
				fprintf(stderr, "Invalid value for sample argument\n");
				return -1;
			}
		} else if (c == 't' || (opt_idx && !strcmp("seed", long_options[opt_idx].name))) {
			seed = atoi(optarg);
		} else if (c == 'l' || (opt_idx && !strcmp("scaleL", long_options[opt_idx].name))) {
			scale = 1;
		} else if (c == 'x' || (opt_idx && !strcmp("approx", long_options[opt_idx].name))) {
			approximate = 1;
		} else if (c == 'h' || (opt_idx && !strcmp("help", long_options[opt_idx].name)))
			printf(HELP_MESSAGE);
		else if (c == 'r' || (opt_idx && !strcmp("restore", long_options[opt_idx].name)))
			restartfile = optarg; // it's fine not to copy, since optarg points into argv
		else {
			printf("unknown option: %c\n\n", c);
			printf(HELP_MESSAGE);
			return -1;
		}
			
	}
	
	if (optind + 2 >= argc) {
		printf(HELP_MESSAGE);
		return -1;
	}
	inputfile = argv[optind];
	dataset = argv[optind+1];
	outputfile = argv[optind+2];

	if (verbose)
		printf("PRNG seed: %d\n", seed);
	srand(seed);
	
	float eps = 1e-3, spl = 0, spz=0.5, lap=1.0;
	int l = 0;
	int n = 0;
	float* tX = readMatrixFromHdf5(inputfile, dataset, &l, &n);
	tX = sample(tX, &l, n, sampleFraction);
	
	//FABIA expects X in column-major format with samples in columns 
	float *tmp = toColumnMajor(tX, l, n);
	free(tX);
	float *X = transposeMatrixCM(tmp, l, n);
	free(tmp);
	
		
	if (verbose)
		printf("data dimension after sampling: (%d, %d)\n", l, n);

	
	// create random initial data
	float* Psi = (float*) malloc(n*sizeof(float));
	float* L = (float*) malloc(n*k*sizeof(float));
	float* Z = (float*) malloc(k*l*sizeof(float));
	float* lapla = (float*) malloc(l*k*sizeof(float));
	
	if (!restartfile) {
		for (int i = 0; i < n*k; ++i) L[i] = (float) rand_normal();
		for (int i = 0; i < l*k; ++i) lapla[i] = 1.0;
		for (int i = 0; i < n; ++i) Psi[i] = 0.2;
		memset(Z, 0, k*l*sizeof(float));
	} else {
		printf("reading startup data from %s\n", restartfile);
		Psi = readMatrixFromRestartFile(restartfile, "Psi", n, 1);
		L = readMatrixFromRestartFile(restartfile, "L", n, k);
		
		// lapla/Z are datapoint specific, and due to sampling we now have different datapoints
		if (sampleFraction < 1.0) {
			memset(Z, 0, k*l*sizeof(float));
			for (int i = 0; i < l*k; ++i) lapla[i] = 1.0;
		} else {
			Z = readMatrixFromRestartFile(restartfile, "Z", k, l);
			lapla = readMatrixFromRestartFile(restartfile, "lapla", l, k);
		}
		if (!Psi || !L || !Z || !lapla)
			return -1;
	}

	if (approximate)
		approx_fabia_cm_f(k, n, l, X, Psi, L, Z, lapla, cyc, (float)alpha, eps, spl, spz, scale, lap, verbose, nthreads);
	else
		fabia_cm_f(k, n, l, X, Psi, L, Z, lapla, cyc, (float)alpha, eps, spl, spz, scale, lap, verbose, nthreads);


	storeResults(outputfile, k, n, l, L, Z, Psi, lapla);
	
	tX = toRowMajor(X, n, l);
	storeMatrixInHdf5(outputfile, "X", tX, n, l);
	free(tX);
	

	
	free(X);
	free(Psi);
	free(L);
	free(Z);
	free(lapla);
	
	return 0;
}
