CC            = gcc
CFLAGS        = -O3 -DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include
CLINKER       = -lgfortran -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread -lm

CXX           = g++
CXXFLAGS      = -O3 -DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include
CXXLINKER     = -lgfortran -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread -lm

NCC           = nvcc
NCFLAGS       = -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -Xcompiler "-DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include" -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse
NCLINKER      = -Xlinker "-lgfortran -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm"

FC            = gfortran
FFLAGS        = -O3 -fdefault-integer-8 -fopenmp -m64 -I${MKLROOT}/include
FLINKER       = -lc -lstdc++ -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread -lm -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse

