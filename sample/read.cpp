#include <iostream>
#include <fstream>
#include <numeric>      // std::partial_sum
#include <complex>
using namespace std;

extern"C" {
void main_cuda_( int *nspin, int *nTerm, double *coeff_lst, int *nbody_lst, int *pos_i_idx, int *pos_i_lst, int *dim_i_lst, int *mat_i_idx, std::complex<double> *mat_i_lst, int *ham_dim, int *nspin_dim, std::complex<double> *v, int *pos_idx, int *mat_idx, int *k, int *maxThreadsPerBlock, int *maxGridSize, int *tn, double *ta, int *m, double *tol, int *itrace, std::complex<double> *w_seq, int *w_seq_len );
}

int main () 
{
    int i;

    int w_seq_len;
    complex<double> * w_seq;

  
    int  klim  = 10;               //  Lanczos factorization length;
    int maxThreadsPerBlock  = 256; //  CUDA grid limitation;
    int maxGridSize[3]      = {2147483647,65535,65535};
    int m                   = 30;  //  zgexpv coefficients;
    double  tol             = 1.0e-12;
    int  itrace          = 0;


    ifstream file ("DipolarDynamics_mat.dat", ios::in|ios::binary);

    int nSpin, nTerm, total_nbody, total_dim, nDim;
    file.read ((char *)&nSpin, sizeof(int) );
    file.read ((char *)&nTerm, sizeof(int) );
    file.read ((char *)&nDim,  sizeof(int) );
    cout << nSpin << "\t" << nTerm << "\t" << nDim <<endl;

    int    * spin_dim   = new int    [nSpin];
    int    * nBody_list = new int    [nTerm];
    double * coeff_list = new double [nTerm];
    int    * pos_offset = new int    [nTerm+1];
    file.read((char *) spin_dim,   nSpin*sizeof(int));
    file.read((char *) coeff_list, nTerm*sizeof(double));
    file.read((char *) nBody_list, nTerm*sizeof(int));

    pos_offset[0]=0;
    partial_sum (nBody_list, nBody_list+nTerm, pos_offset+1);
    total_nbody=pos_offset[nTerm];
    cout << total_nbody << endl;

    int * pos_list   = new int [total_nbody];
    int * dim_list   = new int [total_nbody];
    int * dim2       = new int [total_nbody];
    int * mat_offset = new int [total_nbody+1];
    file.read((char *) pos_list, total_nbody*sizeof(int));
    file.read((char *) dim_list, total_nbody*sizeof(int));
    cout << pos_list[0] << "\t" << pos_list[1] << endl;
    cout << dim_list[0] << "\t" << dim_list[1] << endl;

    for(i=0;i<total_nbody;i++)
        dim2[i]=dim_list[i]*dim_list[i];
    mat_offset[0]=0;
    partial_sum(dim2, dim2+total_nbody, mat_offset+1);
    delete [] dim2;

    total_dim=mat_offset[total_nbody];
    complex<double> * matC;
    double * matR = new double [total_dim];
    double * matI = new double [total_dim];
    matC = new complex<double> [total_dim];
    file.read((char *) matR, total_dim*sizeof(double));
    file.read((char *) matI, total_dim*sizeof(double));
    for(i=0; i<total_dim; i++)
        matC[i] = complex<double> ( matR[i], matI[i] );

    cout << "total_dim=" << total_dim << endl;
    cout << matR[0] << "\t" << matR[1] << "\t" << matR[2] << "\t" << matR[3] << "\t" << matR[4] << endl;
    cout << matI[0] << "\t" << matI[1] << "\t" << matI[2] << "\t" << matI[3] << "\t" << matI[4] << endl;

    file.close();

    double * vecR = new double [nDim];
    double * vecI = new double [nDim];
    ifstream fileVect ("DipolarDynamics_vec.dat", ios::in|ios::binary);
    fileVect.read((char *) vecR, nDim*sizeof(double));
    fileVect.read((char *) vecI, nDim*sizeof(double));
    fileVect.close();

    complex<double> * vecC;
    vecC = new complex<double> [ nDim ];
    for ( i = 0; i < nDim; i++ )
        vecC[i] = complex<double>( vecR[i], vecI[i] );

    int nt;
    double * tlist;
    ifstream fileTime ("DipolarDynamics_time.dat", ios::in|ios::binary);
    fileTime.read((char *) &nt, sizeof(int));

    tlist = new double [nt];
    fileTime.read((char *) tlist, nt*sizeof(double));
    fileTime.close();

    cout << nt << "\t" << tlist[0] << "\t" << tlist[1] << "\t" << tlist[2] << endl;

///////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////

    w_seq_len = nDim * nt;
    w_seq = new complex<double> [w_seq_len];

    
    // expokit in FORTRAN;
    main_cuda_( &nSpin, &nTerm, coeff_list, nBody_list, pos_offset, pos_list, dim_list, mat_offset, matC, &nDim, spin_dim, vecC, &total_nbody, &total_dim, &klim, &maxThreadsPerBlock, maxGridSize, &nt, tlist, &m, &tol, &itrace, w_seq, &w_seq_len );

    ofstream result("res_file.dat", ios::binary);
    double * w_seq_dbl = new double [w_seq_len];

    for ( i = 0; i < w_seq_len; i++ )
        w_seq_dbl[i] = w_seq[i].real();
    result.write((char *)w_seq_dbl , w_seq_len * sizeof(double));

    for ( i = 0; i < w_seq_len; i++ )
        w_seq_dbl[i] = w_seq[i].imag();
    result.write((char *)w_seq_dbl , w_seq_len * sizeof(double));

    result.close();
}
