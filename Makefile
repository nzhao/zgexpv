########################################################################

CUDAAPI_DIR   := /usr/local/cuda/src
ROOT_DIR      := $(shell pwd)

SRC_DIR       := $(ROOT_DIR)/src
OBJ_DIR       := $(ROOT_DIR)/obj
BIN_DIR       := $(ROOT_DIR)/bin

INC_DIR       := $(SRC_DIR)/include

EXPOKIT_DIR   := $(SRC_DIR)/expokit
GPU_DIR       := $(SRC_DIR)/gpu
CPU_DIR       := $(SRC_DIR)/cpu
ZKMV_DIR      := $(SRC_DIR)/zkmv

GPUMAIN_SRC   := evolution_gpu.cpp
CPUMAIN_SRC   := evolution_cpu.cpp

PARAMETER_SRC := para.cpp

EXPOKIT_SRC   := expokit.f mataid.f DLARAN.f
GPU_C_SRC     := hamvec_cuda3.cu
GPU_F_SRC     := hamvec_zgexpv_w_cuda_profile.f main_cuda.f
CPU_C_SRC     := hamvec_func3.cpp
CPU_F_SRC     := hamvec_zgexpv_w_mkl_profile.f main_mkl.f

CUDAAPI_SRC   := fortran.c cusparse_fortran.c

########################################################################

NCC           := nvcc
NCFLAGS       := -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -Xcompiler "-DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include" -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse
NCLINKER      := -Xlinker "-lgfortran -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm"

FC            := gfortran
FFLAGS        := -O3 -fdefault-integer-8 -fopenmp -m64 -I${MKLROOT}/include
FLINKER       := -lc -lstdc++ -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread -lm -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse

########################################################################

SRCDIR_LIST   := $(CUDAAPI_DIR) $(EXPOKIT_DIR) $(SRC_DIR) $(GPU_DIR) $(CPU_DIR)

vpath %.h   $(INC_DIR)
vpath %.c   $(SRCDIR_LIST)
vpath %.cpp $(SRCDIR_LIST)
vpath %.cu  $(SRCDIR_LIST)
vpath %.f   $(SRCDIR_LIST)
vpath %.f90 $(SRCDIR_LIST)

GPUMAIN_OBJ   := $(GPUMAIN_SRC:%.cpp=$(OBJ_DIR)/%.o)
CPUMAIN_OBJ   := $(CPUMAIN_SRC:%.cpp=$(OBJ_DIR)/%.o)
PARAMETER_OBJ := $(PARAMETER_SRC:%.cpp=$(OBJ_DIR)/%.o)
CUDAAPI_OBJ   := $(CUDAAPI_SRC:%.c=$(OBJ_DIR)/%.o)
EXPOKIT_OBJ   := $(EXPOKIT_SRC:%.f=$(OBJ_DIR)/%.o)
GPU_C_OBJ     := $(GPU_C_SRC:%.cu=$(OBJ_DIR)/%.o)
GPU_F_OBJ     := $(GPU_F_SRC:%.f=$(OBJ_DIR)/%.o)
CPU_C_OBJ     := $(CPU_C_SRC:%.cpp=$(OBJ_DIR)/%.o)
CPU_F_OBJ     := $(CPU_F_SRC:%.f=$(OBJ_DIR)/%.o)

GPU_OBJ       := $(CUDAAPI_OBJ) $(EXPOKIT_OBJ) $(GPU_C_OBJ) $(GPU_F_OBJ) $(GPUMAIN_OBJ) $(PARAMETER_OBJ)

CPU_OBJ       := $(CUDAAPI_OBJ) $(EXPOKIT_OBJ) $(CPU_C_OBJ) $(CPU_F_OBJ) $(CPUMAIN_OBJ) $(PARAMETER_OBJ)

all: gpu cpu

gpu: $(GPU_OBJ) | $(BIN_DIR)
	@$(NCC) $(NCFLAGS) $(GPU_OBJ) $(NCLINKER) -o $(BIN_DIR)/gpu_evolve

cpu: $(CPU_OBJ) | $(BIN_DIR)
	@$(NCC) $(NCFLAGS) $(CPU_OBJ) $(NCLINKER) -o $(BIN_DIR)/cpu_evolve

$(GPU_OBJ): | $(OBJ_DIR)

$(CPU_OBJ): | $(OBJ_DIR)

$(GPUMAIN_OBJ):$(OBJ_DIR)/%.o:%.cpp para.h
	@$(NCC) $(NCFLAGS) -I$(INC_DIR) $< -c -o $@

$(CPUMAIN_OBJ):$(OBJ_DIR)/%.o:%.cpp para.h
	@$(NCC) $(NCFLAGS) -I$(INC_DIR) $< -c -o $@

$(PARAMETER_OBJ):$(OBJ_DIR)/%.o:%.cpp para.h
	@$(NCC) $(NCFLAGS) -I$(INC_DIR) $< -c -o $@

para.h: using_namespace.h

$(CUDAAPI_OBJ):$(OBJ_DIR)/%.o:%.c
	@$(NCC) $(NCFLAGS) $< -c -o $@

$(EXPOKIT_OBJ):$(OBJ_DIR)/%.o:%.f
	@$(FC) $(FFLAGS) $< -c -o $@

$(GPU_C_OBJ):$(OBJ_DIR)/%.o:%.cu
	@$(NCC) $(NCFLAGS) $< -c -o $@

$(GPU_F_OBJ):$(OBJ_DIR)/%.o:%.f
	@$(FC) $(FFLAGS) $< -c -o $@

$(CPU_C_OBJ):$(OBJ_DIR)/%.o:%.cpp
	@$(NCC) $(NCFLAGS) $< -c -o $@

$(CPU_F_OBJ):$(OBJ_DIR)/%.o:%.f
	@$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR):
	@mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	@mkdir -p $(BIN_DIR)

.PHONY: clean clean-gpu clean-cpu

clean:
	@rm -rf $(OBJ_DIR)
	@rm -rf $(BIN_DIR)

clean-gpu:
	@rm -f $(GPU_OBJ)
	@rm -f $(BIN_DIR)/gpu.x

clean-cpu:
	@rm -f $(CPU_OBJ)
	@rm -f $(BIN_DIR)/cpu.x

