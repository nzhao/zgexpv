########################################################################

NCC           := nvcc
NCFLAGS       := -arch=sm_30 -O3 -DCUBLAS_GFORTRAN -Xcompiler "-DMKL_ILP64 -fopenmp -m64 -I${MKLROOT}/include" -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse
NCLINKER      := -Xlinker "-lgfortran -L${MKLROOT}/lib/intel64 -lmkl_intel_ilp64 -lmkl_core -lmkl_gnu_thread -ldl -lpthread -lm"

FC            := gfortran
FFLAGS        := -O3 -fdefault-integer-8 -fopenmp -m64 -I${MKLROOT}/include
FLINKER       := -lc -lstdc++ -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_gf_ilp64.a ${MKLROOT}/lib/intel64/libmkl_core.a ${MKLROOT}/lib/intel64/libmkl_gnu_thread.a -Wl,--end-group -ldl -lpthread -lm -L/usr/local/cuda/lib64 -lcudart -lcublas -lcusparse

########################################################################

CUDA_DIR      := /usr/local/cuda
ROOT_DIR      := $(shell pwd)

SRC_DIR       := $(ROOT_DIR)/src
OBJ_DIR       := $(ROOT_DIR)/obj
BIN_DIR       := $(ROOT_DIR)/bin

CUDAAPI_SRC   := fortran.c cusparse_fortran.c

EXPOKIT_SRC   := expokit.f mataid.f DLARAN.f
GPU_C_SRC     := hamvec_cuda3.cu
GPU_F_SRC     := hamvec_zgexpv_w_cuda_profile.f main_cuda.f
CPU_C_SRC     := hamvec_func3.cpp
CPU_F_SRC     := hamvec_zgexpv_w_mkl_profile.f main_mkl.f

EXPOKIT_DIR   := $(SRC_DIR)/expokit
GPU_DIR       := $(SRC_DIR)/gpu
CPU_DIR       := $(SRC_DIR)/cpu
ZKMV_DIR      := $(SRC_DIR)/zkmv

########################################################################

CUDAAPI_LIST  := $(addprefix $(CUDA_DIR)/src/,$(CUDAAPI_SRC))
EXPOKIT_LIST  := $(addprefix $(SRC_DIR)/expokit/,$(EXPOKIT_SRC))
GPU_C_LIST    := $(addprefix $(SRC_DIR)/gpu/,$(GPU_C_SRC))
GPU_F_LIST    := $(addprefix $(SRC_DIR)/gpu/,$(GPU_F_SRC))
CPU_C_LIST    := $(addprefix $(SRC_DIR)/cpu/,$(CPU_C_SRC))
CPU_F_LIST    := $(addprefix $(SRC_DIR)/cpu/,$(CPU_F_SRC))

CUDAAPI_OBJ   := $(addprefix $(OBJ_DIR)/,$(patsubst %.c,%.o,$(CUDAAPI_SRC)))
EXPOKIT_OBJ   := $(addprefix $(OBJ_DIR)/,$(patsubst %.f,%.o,$(EXPOKIT_SRC)))
GPU_C_OBJ     := $(addprefix $(OBJ_DIR)/,$(patsubst %.cu,%.o,$(GPU_C_SRC)))
GPU_F_OBJ     := $(addprefix $(OBJ_DIR)/,$(patsubst %.f,%.o,$(GPU_F_SRC)))
CPU_C_OBJ     := $(addprefix $(OBJ_DIR)/,$(patsubst %.cpp,%.o,$(CPU_C_SRC)))
CPU_F_OBJ     := $(addprefix $(OBJ_DIR)/,$(patsubst %.f,%.o,$(CPU_F_SRC)))

GPU_OBJ       := $(CUDAAPI_OBJ) $(EXPOKIT_OBJ) $(GPU_C_OBJ) $(GPU_F_OBJ) $(OBJ_DIR)/evolution_gpu.o $(OBJ_DIR)/para.o

CPU_OBJ       := $(CUDAAPI_OBJ) $(EXPOKIT_OBJ) $(CPU_C_OBJ) $(CPU_F_OBJ) $(OBJ_DIR)/evolution_cpu.o $(OBJ_DIR)/para.o

echo: $(OBJ_DIR)/cusparse_fortran.o
	@echo $(EXPOKIT_SRC)
	@echo $(EXPOKIT_LIST)
	@echo $(EXPOKIT_OBJ)
	@echo $(CUDAAPI_SRC)
	@echo $(CUDAAPI_LIST)
	@echo $(CUDAAPI_OBJ)
	@echo $(GPU_OBJ)

gpu: $(GPU_OBJ) | $(BIN_DIR)
	@$(NCC) $(NCFLAGS) $(GPU_OBJ) $(NCLINKER) -o $(BIN_DIR)/gpu.x

cpu: $(CPU_OBJ) | $(BIN_DIR)
	@$(NCC) $(NCFLAGS) $(CPU_OBJ) $(NCLINKER) -o $(BIN_DIR)/cpu.x

$(GPU_OBJ): | $(OBJ_DIR)

$(CPU_OBJ): | $(OBJ_DIR)

$(OBJ_DIR)/cusparse_fortran.o: $(CUDA_DIR)/src/cusparse_fortran.c
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR)/fortran.o: $(CUDA_DIR)/src/fortran.c
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR)/expokit.o: $(EXPOKIT_DIR)/expokit.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/mataid.o: $(EXPOKIT_DIR)/mataid.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/DLARAN.o: $(EXPOKIT_DIR)/DLARAN.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/hamvec_cuda3.o: $(GPU_DIR)/hamvec_cuda3.cu
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR)/hamvec_zgexpv_w_cuda_profile.o: $(GPU_DIR)/hamvec_zgexpv_w_cuda_profile.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/main_cuda.o: $(GPU_DIR)/main_cuda.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/evolution_gpu.o: $(SRC_DIR)/evolution_gpu.cpp
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR)/para.o: $(SRC_DIR)/para.cpp
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR)/hamvec_func3.o: $(CPU_DIR)/hamvec_func3.cpp
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR)/hamvec_zgexpv_w_mkl_profile.o: $(CPU_DIR)/hamvec_zgexpv_w_mkl_profile.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/main_mkl.o: $(CPU_DIR)/main_mkl.f
	$(FC) $(FFLAGS) $< -c -o $@

$(OBJ_DIR)/evolution_cpu.o: $(SRC_DIR)/evolution_cpu.cpp
	$(NCC) $(NCFLAGS) $< -c -o $@

$(OBJ_DIR):
	mkdir -p $(OBJ_DIR)

$(BIN_DIR):
	mkdir -p $(BIN_DIR)

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

