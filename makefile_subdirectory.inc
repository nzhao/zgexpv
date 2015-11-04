include $(ROOT_DIR)/makefile_compiler.inc

########################################################################
##Makefile for subdirectory
########################################################################

SUBDIRS       = $(shell ls -l | grep ^d | awk '{print $$9}')

SOURCE_C      = ${wildcard *.c}
OBJS_C        = ${patsubst %.c, %.o, $(SOURCE_C)}

SOURCE_CPP    = ${wildcard *.cpp}
OBJS_CPP      = ${patsubst %.cpp, %.o, $(SOURCE_CPP)}

SOURCE_CU     = ${wildcard *.cu}
OBJS_CU       = ${patsubst %.cu, %.o, $(SOURCE_CU)}

SOURCE_F      = ${wildcard *.f}
OBJS_F        = ${patsubst %.f, %.o, $(SOURCE_F)}

SOURCE_F90    = ${wildcard *.f90}
OBJS_F90      = ${patsubst %.f90, %.o, $(SOURCE_F90)}

all: $(SUBDIRS) $(OBJS_C) $(OBJS_CPP) $(OBJS_CU) $(OBJS_F) $(OBJS_F90)

$(SUBDIRS): ECHO
	@make -C $@

$(OBJS_C): %.o: %.c
	@$(CC) $(CFLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@

$(OBJS_CPP): %.o: %.cpp
	@$(CXX) $(CXXFLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@

$(OBJS_CU): %.o: %.cu
	@$(NCC) $(NCFLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@

$(OBJS_F): %.o: %.f
	@$(FC) $(FFLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@

$(OBJS_F90): %.o: %.f90
	@$(FC) $(FFLAGS) -c $^ -o $(ROOT_DIR)/$(OBJS_DIR)/$@

ECHO:
	@echo $(SUBDIRS)
########################################################################

