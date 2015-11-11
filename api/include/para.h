#ifndef PARA_H_
#define PARA_H_

#include "using_namespace.h"

extern char INPUT_PATH[500];
extern char INPUT_FILE[500];
extern char OUTPUT_PATH[500];
extern char OUTPUT_FILE[500];

extern double KRYLOV_TOL;
extern unsigned long KRYLOV_M;
extern unsigned long KRYLOV_ITRACE;
extern unsigned long KLIM;
extern unsigned long MAX_THREADS_PER_BLOCK;
extern unsigned long MAX_GRID_SIZE_1;
extern unsigned long MAX_GRID_SIZE_2;
extern unsigned long MAX_GRID_SIZE_3;

void Load_Config();
void ParameterResolve(int argc, char ** argv);
void PrintParameters();

#endif /* PARAMETERS_H_ */
