#This script compiles the fortran code and generates a python module
export PATH=/share/anaconda3/bin/:$PATH
export LD_LIBRARY_PATH=/share/anaconda3/lib/:$LD_LIBRARY_PATH

export FC=gfortran
export F90=gfortran

#Optimized
f2py  -c -lgomp --f90flags="-fopenmp -lgomp -O3" -m common_functions common_functions.f90

#Debug
#f2py  -c --f90flags="-g -traceback" -m common_functions common_functions.f90
