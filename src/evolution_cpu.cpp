#include <iostream>
#include <fstream>
#include <numeric>      // std::partial_sum
#include <complex>
#include "include/para.h"
using namespace std;

int engine();
extern"C" {
void main_mkl_( size_t *nspin, size_t *nTerm, double *coeff_lst, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *mat_i_lst, size_t *ham_dim, size_t *nspin_dim, std::complex<double> *v, size_t *pos_idx, size_t *mat_idx, size_t *k, size_t *tn, double *ta, size_t *m, double *tol, size_t *itrace, std::complex<double> *w_seq, size_t *w_seq_len );
}

int main (int argc, char **argv) 
{

    ParameterResolve(argc, argv);
    PrintParameters();
    engine();
}

int engine()
{
    size_t i;
  
    size_t m = KRYLOV_M;
    double tol = KRYLOV_TOL;
    size_t itrace = KRYLOV_ITRACE;
    size_t klim = KLIM;
    size_t maxThreadsPerBlock = MAX_THREADS_PER_BLOCK;
    size_t maxGridSize[3];
    maxGridSize[0] = MAX_GRID_SIZE_1;
    maxGridSize[1] = MAX_GRID_SIZE_2;
    maxGridSize[2] = MAX_GRID_SIZE_3;

    size_t  nSpin, nTerm, nDim, total_nbody, total_dim;
    size_t  *nBody_list, *pos_offset, *pos_list, *dim_list, *mat_offset, *spin_dim;
    double * coeff_list;
    

//////////////////////////////////////////////////////////////////////////
//read Matrix
    char filename_in[500];
    strcpy(filename_in, INPUT_PATH);
    strcat(filename_in, INPUT_FILE);
    ifstream file (filename_in, ios::in|ios::binary);
    
    total_nbody = 0; total_dim = 0; nSpin = 0; nTerm = 0; nDim = 0;
    file.read ((char *)&nSpin, sizeof(size_t) );
    file.read ((char *)&nTerm, sizeof(size_t) );
    file.read ((char *)&nDim,  sizeof(size_t) );
    
    spin_dim   = new size_t [nSpin];
    nBody_list = new size_t [nTerm];
    coeff_list = new double [nTerm];
    pos_offset = new size_t [nTerm+1];
    file.read((char *) spin_dim,   nSpin*sizeof(size_t));
    file.read((char *) coeff_list, nTerm*sizeof(double));
    file.read((char *) nBody_list, nTerm*sizeof(size_t));

    pos_offset[0]=0;
    partial_sum (nBody_list, nBody_list+nTerm, pos_offset+1);
    total_nbody=pos_offset[nTerm];

    pos_list   = new size_t [total_nbody];
    dim_list   = new size_t [total_nbody];
    mat_offset = new size_t [total_nbody+1];
    file.read((char *) pos_list, total_nbody*sizeof(size_t));
    file.read((char *) dim_list, total_nbody*sizeof(size_t));
    
    size_t * dim2 = new size_t [total_nbody];
    for(i=0;i<total_nbody;i++)
        dim2[i]=dim_list[i]*dim_list[i];
    mat_offset[0]=0;
    partial_sum(dim2, dim2+total_nbody, mat_offset+1);
    delete [] dim2;
    total_dim=mat_offset[total_nbody];
    
    complex<double> * matC;
    double * matR, * matI;

    matR = new double [total_dim];
    matI = new double [total_dim];
    matC = new complex<double> [total_dim];
    file.read((char *) matR, total_dim*sizeof(double));
    file.read((char *) matI, total_dim*sizeof(double));
    for(i=0; i<total_dim; i++)
        matC[i] = complex<double> ( matR[i], matI[i] );

    //////////////////////////////////////////////////////////////////////////
    //read vector

    double * vecR, * vecI;
    complex<double> * vecC;

    vecR = new double [nDim];
    vecI = new double [nDim];
    vecC = new complex<double> [ nDim ];

    file.read((char *) vecR, nDim*sizeof(double));
    file.read((char *) vecI, nDim*sizeof(double));

    for ( i = 0; i < nDim; i++ )
        vecC[i] = complex<double>( vecR[i], vecI[i] );

    //////////////////////////////////////////////////////////////////////////
    //read time seq

    size_t nt=0;
    double * tlist;
    file.read((char *) &nt, sizeof(size_t));

    tlist = new double [nt];
    file.read((char *) tlist, nt*sizeof(double));

///////////////////////////////////////////////////////////////////////////////////////////
// Expv computation

    size_t w_seq_len;
    complex<double> * w_seq;
    
    w_seq_len = nDim * nt;
    w_seq = new complex<double> [w_seq_len];

    // expokit in FORTRAN;
    main_mkl_( &nSpin, &nTerm, coeff_list, nBody_list, pos_offset, pos_list, dim_list, mat_offset, matC, &nDim, spin_dim, vecC, &total_nbody, &total_dim, &klim, &nt, tlist, &m, &tol, &itrace, w_seq, &w_seq_len );

///////////////////////////////////////////////////////////////////////////////////////////
// output result

    char filename_out[500];
    strcpy(filename_out, OUTPUT_PATH);
    strcat(filename_out, OUTPUT_FILE);
    ofstream result(filename_out, ios::binary);

    double * w_seq_dbl = new double [w_seq_len];

    cout << nDim << "exported." << endl;
    result.write((char *)&nDim, sizeof(size_t));
    result.write((char *)&nt,   sizeof(size_t));
    cout << nt  << "exported." << endl;

    for ( i = 0; i < w_seq_len; i++ )
        w_seq_dbl[i] = w_seq[i].real();
    result.write((char *)w_seq_dbl , w_seq_len * sizeof(double));

    for ( i = 0; i < w_seq_len; i++ )
        w_seq_dbl[i] = w_seq[i].imag();
    result.write((char *)w_seq_dbl , w_seq_len * sizeof(double));

    result.close();
    delete [] w_seq_dbl;
    
    return 0;
}
