
##===========================================================

## Minimal Makefile for Fortran parts of KronLinInv


##
## Note that Julia and Python have their own "pure" implementation of the algorithm,
##  which doesn't need any compilation. The respective Fortran wrappers are optional,
##  for speed.
##
## It is possible to compile single versions of the program or all together:
##
## * sharedmemfortran:  compile the serial/shared memory Fortran code
## * distribmemfortran: compile the OpenMPI/distributed memory Fortran code
## * wrapforjulia:      compile Fortran wrappers for Julia
## * wrapforpython:     compile Fortran wrappers for Python
## * all:               compile all Fortran code and wrappers for Julia and Python
##
##  Example:
##  make distribmemfortran
##  will compile the OpenMPI distributed memory version of KronLinInv		
##  make all
##  will compile all Fortran files
##



## MPI Fortran compiler
MPIFC = mpifort
## Fortran compiler
FC = gfortran

## Output dir for Fortran executables
OUTDIR = examples/

## Shared memory parallelism, set to nothing if OpenMP is not available 
OMP = -fopenmp 

## options for gfortran, may differ for other compilers
OPTFOR =  -std=f2008 -O3 -funroll-loops -llapack -lblas

## DEBUG
##OPTFOR = -std=f2008 -g -fcheck=all -fbounds-check -Wrealloc-lhs-all -llapack -lblas

## include/linking for HDF5 libraries, system-depent
LIBSHDF = -I/usr/include/hdf5/serial -L/usr/lib/x86_64-linux-gnu/hdf5/serial



##===========================================================

all: sharedmemfortran distribmemfortran wrapforjulia wrapforpython



## Compile test program in Fortran, serial/shared memory
sharedmemfortran:
	$(FC)  -std=f2008 -O3 -fno-realloc-lhs fortran/kronlininv.f08 fortran/rdwrhdf5.f08  fortran/test_kronlininv.f08 $(LIBSHDF) -lblas -llapack -lhdf5 -lhdf5_fortran $(OMP) -o $(OUTDIR)test.x 
	rm *.mod #*.o *.so # clean module and obj files



## Compile test program in Fortran, distributed memory (OpenMPI)
distribmemfortran: 
	$(MPIFC) -std=f2008 -O3 -fno-realloc-lhs fortran/ompi_kronlininv.f08 fortran/rdwrhdf5.f08  fortran/test_ompi_kronlininv.f08 $(LIBSHDF) -lblas -llapack -lhdf5 -lhdf5_fortran -o $(OUTDIR)testMPI.x 
	rm *.mod  # clean moduue and obj files



## Compile Fortran wrappers for Julia
wrapforjulia: 
	$(FC)  -std=f2008 -O3 -funroll-loops -lblas -llapack fortran/kronlininv.f08  julia/wrap_kronlininv_julia.f90 $(OMP) -fPIC -shared -o julia/kronlininv_wrapjulia.so
	rm *.mod # clean module and obj files



## Compile Fortran wrappers for Python using f2py
wrapforpython: 
# echo "Creating double-type mapping for f2py in file .f2py_f2cmap"
# cat <<EOF > python/.f2py_f2cmap
# dict( real= dict(dp='double', c_double='double'),integer = dict(c_int='int') )
# EOF
	$(FC) -std=f2008 -O3 -fPIC -c fortran/kronlininv.f08 -o python/kronlininv.o $(OMP)

	f2py --f90flags="-fopenmp" -lgomp  -llapack -lblas  -c python/kronlininv.o  python/wrap_kronlininv_f2py.f90  -m kronlininv_f2py

	mv kronlininv_f2py.so python/kronlininv_f2py.so
	rm *.mod python/*.o # clean module and obj files



