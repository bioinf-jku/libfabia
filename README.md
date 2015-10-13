libfabia
========

(c) 2013 Thomas Unterthiner
Last Change: 2013-11-07 (Version 0.03)

This library implements the FABIA biclustering algorithm from


*FABIA: Factor Analysis for Bicluster Acquisition*, Sepp Hochreiter et al,,  Bioinformatics 2010

This implementation aims at HPC machines, parallelizing fabia over as many cores as possible to make it feasible to extract hundreds or thousands of biclusters. The program has been tested with using up to 500 CPU cores in parallel, and the speedups are almost linear.

For Input/Output, this program uses HDF5 files. Meaning it expects its input in an HDF5 file (you'll have to pass both file and dataset name as arguments) and will store its final matrices (L, Z, Psi, Lapla) in HDF5 as well.


Features
--------

- Uses LAPACK/BLAS for its numercial operations.
- Multithreaded with linear speedup.
- Creates rollback points after a user-defined number of iterations.
  Calculations can restart from such rollback-points (or, since rollback
  points are simply intermediary result-files, also from previous runs of
  the algorithm).
- Detects and re-initializes empty bi-clusters that may appear during runs.
- can interface with R


Usage
-----

``fabia inputfile datasetname outputfile``

Parameters:

 -v/--verbose n     print progress status and create rollback point every n cycles.
 -h/--help          show help message.
 -k/--nbicluster k  look for n biclusters.
 -n/--nthreads nt   use a total of nt threads.
 -c/--cyc  c        run FABIA for c cycles.
 -a/--alpha  a      set Laplace prior to a.
 -s/--sample f      sample fraction f from the data, instead of using all of it.
 -r/--restart fn    restart computation using results from file fn.



Compilation instructions
------------------------

The provided Makefile should work out of the box, if all necessary libraries are installed on the system. Typing ``make all`` should do the trick.

The package requires libhdf5, as well as a BLAS and a LAPACK library. Additionally the common fortran-compatibility headers for LAPACK and BLAS (f2c.h clapack.h blas.h) are needed, and should be placed in ext/include (just google them, they're standard headers and easy to find).



License
-------
libfabia is licensed under the [General Public License (GPL) Version 2 or higher](http://www.gnu.org/licenses/gpl-2.0.html) See ``License.txt`` for details.
