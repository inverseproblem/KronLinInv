# KronLinInv

Kronecker-product-based linear inversion of geophysical (or other kinds of) data under Gaussian and separability assumptions. The code computes the posterior mean model and the posterior covariance matrix (or subsets of it) in an efficient manner (parallel algorithm) taking into account 3-D correlations both in the model parameters and in the observed data.

If you use this code for research or else, please cite the related paper:
 
Andrea Zunino, Klaus Mosegaard,
**An efficient method to solve large linearizable inverse problems under Gaussian and separability assumptions**,
Computers & Geosciences, 2018
ISSN 0098-3004, <https://doi.org/10.1016/j.cageo.2018.09.005>.

---

**New**: a pure parallel Julia implementation is available, see <https://github.com/inverseproblem/KronLinInv.jl>

---

Few versions of the code for different languages are provided:

- Fortran (2008): a simple serial version, where OpenMP (parallelisation for shared memory architecture) can be turned on at compile time and a parallel distributed memory OpenMPI version (the "main" implementation);

- Julia: a new package available, see <https://github.com/inverseproblem/KronLinInv.jl>

- Python: a pure Python implementation (quite slow, mostly for learning purposes) and wrappers from Python to the Fortran routines.

The reference implementation of this program is in the form of Fortran subroutines, to be included in the user's code, however, a couple of test programs reading input and writing output from/to HDF5 files are provided for convenience. Please see the documentation attached to the Fortran version at <https://inverseproblem.github.io/KronLinInv/>.



### Prerequisites

- Fortran version:
   - required: LAPACK
   - optional (but recommended): OpenMPI, HDF5

- Julia: Julia version >=1.0

- Python: Python version >=3.?
   - optional: f2py for wrapping Fortran routines

### Installing

Documentation of Julia and Python versions is work in progress.  Please refer to the Doxygen documentation of the Fortran code for now (<https://inverseproblem.github.io/KronLinInv/>). Only the three main subroutines of the Fortran code have been ported to Julia and Python.


#### Fortran

See the Makefile for compiling the Fortran subroutines and the Doxygen documentation <https://inverseproblem.github.io/KronLinInv/>.
The Fortran code is intended to be inserted in the user code in the form of a set of subroutines, however, a couple of test programs calling these subroutines are provided, which read input and write the final output to HDF5 files. The Makefile shows how to compile such programs. 

_The **first** thing to compute is always the set of "factors" using the subroutine calcfactors(), then the posterior mean or covariance (or part of it) can be computed._

_Remark_: The OpenMPI version of the Fortran subroutines needs at least an MPI\_init() before and MPI\_finalize() after from the calling code. A convenience module "parautil" is provided with several utility routines. See "test\_ompi\_kronlininv.f03" for an example, where para\_init() and para\_finish() are used to open and close the MPI interface.

The input HDF5 file for the test programs must contain the following arrays: 
- G1,G2,G3, the 3 forward operators
- Cm1,Cm2,Cm3, the 3 covariance matrices for the model parameters
- Cd1,Cd2,Cd3, the 3 covariance matrices for the observed data
- dobs, the observed data vector
- mprior, the mean prior model vector

The directory "examples" contains an example data set that can be used to run the two Fortran test programs. Compile one or both of the Fortran versions of KronLinInv using the Makefile, place the executable in the "example" directory and then launch the programs with  
./test.x  
or  
mpirun -np [num-of-procs] testMPI.x  
for the shared and distributed memory, respectively.

#### Julia

*What follows is an outdated version of the code, please use the new Julia package KronLinInv.jl, see <https://github.com/inverseproblem/KronLinInv.jl>.*

To use the following functions the file "kronlininv.jl" in the julia directory must be included in your code. For instance: include("kronlininv.jl")
_The **first** thing to compute is always the set of "factors" using the function calcfactors(), then the posterior mean or covariance (or part of it) can be computed._

List of available functions:

- calcfactors(...): Calculate the "factors" needed for subsequent calculations.

- posteriormean(...): Calculate the posterior mean model given input data and the previously calculated "factors" (see calcfactors(...)). Serial and parallel version in pure Julia (call it using posteriormean(..., runparallel=true/false)). 

- blockpostcov(...):  Calculate a block (or all) of the posterior covariance by specifying the indices of the first (astart,aend) and second (bstart,bend) dimensions.

Fortran wrappers (see Makefile for compiling the Fortran code):

- fortran\_calcfactors(...), fortran\_posteriormean(...) and  fortran\_blockpostcov(...): perform the same calculations of the function described above but by calling the Fortran subroutines. Note that fortran\_posteriormean(...) calls the serial Fortran version. The "fortranlibfl" argument shold point to the shared library "kron\_utils.so", for instance: fortranlibfl = "julia/kron\_utils.so".
 

#### Python

To use the following functions import the file kronlininv.py as a module, i.e., "import kronlininv".
_The **first** thing to compute is always the set of "factors" using the function calcfactors(), then the posterior mean or covariance (or part of it) can be computed._

List of available functions:

- calcfactors(...): Calculate the "factors" needed for subsequent calculations.

- posteriormean(...): Calculate the posterior mean model given input data and the previously calculated "factors" (see calcfactors(...)). 

- blockpostcov(...):  Calculate a block (or all) of the posterior covariance by specifying the indices of the first (astart,aend) and second (bstart,bend) dimensions.

Fortran wrappers (see Makefile for compiling the Fortran code):

- fortran\_calcfactors(...), fortran\_posteriormean(...) and  fortran\_blockpostcov(...): perform the same calculations of the function described above but by calling the Fortran subroutines. Note that fortran\_posteriormean(...) calls the serial Fortran version.


## Authors
Andrea Zunino, 
Niels Bohr Institute
