########################################################################
#Makefile for main directory
########################################################################

CUDA_DIR      = /usr/local/cuda

CUDA_OBJS     = fortran.o cusparse_fortran.o

SAMPLE_DIR    = sample

SUBDIRS       = $(shell ls -l | grep ^d | awk '($$9 != "debug")&&($$9 != "data")&&($$9 != "sample") {print $$9}')

ROOT_DIR      = $(shell pwd)

OBJS_DIR      = debug/obj

BIN_DIR       = debug/bin

export CC BIN OBJS_DIR BIN_DIR ROOT_DIR

include $(ROOT_DIR)/makefile_compiler.inc

cuda: all $(SAMPLE_DIR)/main_cuda.cu
	@$(NCC) $(NCFLAGS) $(SAMPLE_DIR)/main_cuda.cu $(OBJS_DIR)/*.o $(NCLINKER) -o $(ROOT_DIR)/$(BIN_DIR)/cuda.x

mkl: all $(SAMPLE_DIR)/main_mkl.cu
	@$(NCC) $(NCFLAGS) $(SAMPLE_DIR)/main_mkl.cu $(OBJS_DIR)/*.o $(NCLINKER) -o $(ROOT_DIR)/$(BIN_DIR)/mkl.x

newCUDA: all $(SAMPLE_DIR)/read.cpp
	@$(NCC) $(NCFLAGS) $(SAMPLE_DIR)/read.cpp $(OBJS_DIR)/*.o $(NCLINKER) -o $(ROOT_DIR)/$(BIN_DIR)/newCUDA.x

all: $(SUBDIRS) $(CUDA_OBJS)

$(SUBDIRS): ECHO
	@make -C $@

ECHO:
	@echo $(SUBDIRS)

$(CUDA_OBJS): %.o: $(CUDA_DIR)/src/%.c
	@$(NCC) $(NCFLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@

.PHONY: clean

clean:
	@rm -f $(OBJS_DIR)/*.o
	@rm -f $(BIN_DIR)/*.x
########################################################################

