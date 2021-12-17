#This script compiles the fortran code and generates a python module
export PATH=$HOME/data/software/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=$HOME/data/software/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=ifort
export F77=ifort
export F90=ifort

f2py  -c -lgomp --f90flags="-qopenmp -lgomp -O3" -m common_functions common_functions.f90

#f2py3  -c --f90flags="-g" -m common_functions common_functions.f90
