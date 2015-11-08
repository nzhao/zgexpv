#include <iostream>
#include <fstream>
#include <numeric>      // std::partial_sum
using namespace std;

int main () 
{
    int i;

    ifstream file ("DipolarDynamics_mat.dat", ios::in|ios::binary);

    int nSpin, nTerm, total_nbody, total_dim, nDim;
    file.read ((char *)&nSpin, sizeof(int) );
    file.read ((char *)&nTerm, sizeof(int) );
    file.read ((char *)&nDim,  sizeof(int) );
    cout << nSpin << "\t" << nTerm <<endl;

    int    * spin_dim   = new int    [nSpin];
    int    * nBody_list = new int    [nTerm];
    double * coeff_list = new double [nTerm+1];
    int    * pos_offset = new int    [nTerm];
    file.read((char *) spin_dim,   nSpin*sizeof(int));
    file.read((char *) nBody_list, nTerm*sizeof(int));
    file.read((char *) coeff_list, nTerm*sizeof(double));

    pos_offset[0]=0;
    std::partial_sum (nBody_list, nBody_list+nTerm, pos_offset+1);
    total_nbody=pos_offset[nTerm];

    int * pos_list   = new int [total_nbody];
    int * dim_list   = new int [total_nbody];
    int * dim2       = new int [total_nbody];
    int * mat_offset = new int [total_nbody+1];
    file.read((char *) pos_list, total_nbody*sizeof(int));
    file.read((char *) dim_list, total_nbody*sizeof(int));

    for(i=0;i<total_nbody;i++)
        dim2[i]=dim_list[i]*dim_list[i];
    mat_offset[0]=0;
    std::partial_sum(dim2, dim2+total_nbody, mat_offset+1);
    delete [] dim2;

    total_dim=mat_offset[total_nbody];
    double * matR = new double [total_dim];
    double * matI = new double [total_dim];
    file.read((char *) matR, total_dim*sizeof(double));
    file.read((char *) matI, total_dim*sizeof(double));


}
