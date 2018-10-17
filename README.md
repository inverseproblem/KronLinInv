# KronLinInv

 Kronecker-product-based linear inversion under Gaussian and separability assumptions.
 
 If you use this code for research or else, please cite the related paper:
 
Andrea Zunino, Klaus Mosegaard,
**An efficient method to solve large linearizable inverse problems under Gaussian and separability assumptions**,
Computers & Geosciences, 2018
ISSN 0098-3004, <https://doi.org/10.1016/j.cageo.2018.09.005>.


Few versions of the code for different languages are provided:

- Fortran (2008): a simple serial version, where OpenMP (parallelisation for shared memory architecture) can be turned on at compile time and a parallel distributed memory OpenMPI version (the "main" implementation);

- Julia: a pure Julia implementation (serial and parallel - fairly fast) and wrappers from Julia to the Fortran routines;

- Python: a pure Python implementation (quite slow, mostly for learning purposes) and wrappers from Python to the Fortran routines.

The reference implementation of this program is in the form of Fortran subroutines, to be included in the user's code, however, a couple of test programs reading input and writing output from/to HDF5 files are provided for convenience. Please see the documentation attached to the Fortran version at https://inverseproblem.github.io/KronLinInv/.



### Prerequisites

- Fortran version:
   - required: LAPACK
   - optional (but recommended): OpenMPI, HDF5

- Julia: Julia version >=1.0

- Python: Python version >=3.?


### Installing


to be written...


## Authors
Andrea Zunino, 
Niels Bohr Institute
