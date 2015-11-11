#include <iostream>
#include <fstream>
#include <complex>
#include <string>
#include <cmath>

#include "../src/ham.cpp"

using namespace std;

void cuzkmv( int nspin, int nTerm, double *coeff_lst, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *mat_i_lst, size_t vlen, size_t *nspin_dim, std::complex<double> *v, size_t pos_idx, size_t mat_idx, std::complex<double> *w );

int main()
{
  //  data input;
  char    ham_input_name[]    = "../../data/input_ham.dat";
  char    vec_input_name[]    = "../../data/input_vec.dat";
  char    tim_input_name[]    = "../../data/input_tim.dat";
  size_t  numbering_base      = 1;
  
  //  data output;
  char    output_name[]       = "../../data/output_cuda.dat";
  std::complex<double> *w;
  
//======================================================================
// hamvec_read
//======================================================================
  size_t  nspin;
  size_t  nTerm;
  double  coeff;
  size_t  nbody;
  size_t  pos_i;
  size_t  dim_i;
  std::complex<double> *mat_i;
  std::complex<double> *v;
  
  size_t  d0, d2;
  
  QuanObj quanobj;
  HamObj  hamobj;
  Ham     ham;
  
  double          *coeff_lst;
  size_t          *nbody_lst;
  size_t          *pos_i_idx, pos_idx;
  size_t          *pos_i_lst;
  size_t          *dim_i_lst;
  size_t          *mat_i_idx, mat_idx;
  std::complex<double> *mat_i_lst;
  size_t          *nspin_dim;
  size_t          ham_dim;
  
  FILE    *fp;
  size_t  i, j, k;
  double  *mat_i_real, *mat_i_imag;
  
  // input hamiltonian;
  fp = fopen( ham_input_name, "rb" );
  if ( fread( &nspin, sizeof(nspin), 1, fp ) != 1 )
    std::cout << "file read error: nspin" << std::endl;
  if ( fread( &nTerm, sizeof(nTerm), 1, fp ) != 1 )
    std::cout << "file read error: nTerm" << std::endl;
  ham.set_nspin( nspin );
  d0 = 0;
  for ( i = 0; i < nTerm; i++ )
  {
    hamobj = HamObj();
    if ( fread( &coeff, sizeof(coeff), 1, fp ) != 1 )
      std::cout << "file read error: coeff" << std::endl;
    hamobj.set_coe( coeff );
    if ( fread( &nbody, sizeof(nbody), 1, fp ) != 1 )
      std::cout << "file read error: nbody" << std::endl;
    for ( j = 0; j < nbody; j++ )
    {
      quanobj = QuanObj();
      if ( fread( &pos_i, sizeof(pos_i), 1, fp ) != 1 )
        std::cout << "file read error: pos_i" << std::endl;
      if ( fread( &dim_i, sizeof(dim_i), 1, fp ) != 1 )
        std::cout << "file read error: dim_i" << std::endl;
      d2 = dim_i * dim_i;
      if ( d0 < d2 )
      {
        d0 = d2;
        mat_i_real = new double [ d2 ];
        mat_i_imag = new double [ d2 ];
        mat_i = new std::complex<double> [ d2 ];
      }
      if ( fread( mat_i_real, sizeof( mat_i_real[0] ), d2, fp ) != d2 )
        std::cout << "file read error: mat_i_real" << std::endl;
      if ( fread( mat_i_imag, sizeof( mat_i_imag[0] ), d2, fp ) != d2 )
        std::cout << "file read error: mat_i_imag" << std::endl;
      for ( k = 0; k < d2; k++ )
        mat_i[k] = std::complex<double> ( mat_i_real[k], mat_i_imag[k] );
      quanobj.set_opr( dim_i, mat_i );
      hamobj.append_obj( 1, &pos_i, &quanobj );
    }
    ham.append_obj( 1, &hamobj );
  }
  fclose( fp );
  //====================================================================
  // regroup data for computation;
  //====================================================================
  // STEP 1: get nTerm -sized lists of [ coeff_lst, nbody_lst and pos_i_idx ], and count for [ pos_idx ] which is the size of pos_i_lst, dim_i_lst, mat_i_idx;
  coeff_lst = new double [ nTerm ];
  nbody_lst = new size_t [ nTerm ];
  pos_i_idx = new size_t [ nTerm ];
  pos_idx = 0;
  for ( i = 0; i < nTerm; i++ )
  {
    hamobj = *ham.get_obj( 1, &i );
    coeff = hamobj.get_coe();
    nbody = hamobj.get_nbody();
    
    coeff_lst[ i ] = coeff;
    nbody_lst[ i ] = nbody;
    pos_i_idx[ i ] = pos_idx;
    
    pos_idx += nbody;
  }
  // STEP 2: get the pos_idx -sized lists of [ pos_i_lst, dim_i_lst, mat_i_idx ], and count for [ mat_idx ] which is the size of mat_i_lst;
  pos_i_lst = new size_t [ pos_idx ];
  dim_i_lst = new size_t [ pos_idx ];
  mat_i_idx = new size_t [ pos_idx ];
  mat_idx = 0;
  for ( i = 0; i < nTerm; i++ )
  {
    hamobj = *ham.get_obj( 1, &i );
    nbody = hamobj.get_nbody();
    for ( j = 0; j < nbody; j++ )
    {
      quanobj = *hamobj.get_obj( 1, &j );
      dim_i = quanobj.get_dim();
      
      pos_i_lst[ pos_i_idx[ i ] + j ] = *hamobj.get_pos( 1, &j ) - numbering_base;
      dim_i_lst[ pos_i_idx[ i ] + j ] = quanobj.get_dim();
      mat_i_idx[ pos_i_idx[ i ] + j ] = mat_idx;
      
      mat_idx += dim_i * dim_i;
    }
  }
  // STEP 3: get the mat_idx -sized list [ mat_i_lst ];
  mat_i_lst = new std::complex<double> [ mat_idx ];
  for ( i = 0; i < nTerm; i++ )
  {
    hamobj = *ham.get_obj( 1, &i );
    nbody = hamobj.get_nbody();
    
    for ( j = 0; j < nbody; j++ )
    {
      quanobj = *hamobj.get_obj( 1, &j );
      dim_i = quanobj.get_dim();
      mat_i = quanobj.get_opr();
      d2 = dim_i * dim_i;
      
      for ( k = 0; k < d2; k++ )
        mat_i_lst[ mat_i_idx[ pos_i_idx[ i ] + j ] + k] = mat_i[ k ];
    }
  }
  //====================================================================
  
  // output the hamiltonian;
  nspin_dim = new size_t [nspin];
  for ( i = 0; i < nspin; i++ )
    nspin_dim[i] = 0;
  for ( i = 0; i < pos_idx; i++ )
  {
    if ( pos_i_lst[i] >= nspin )
    {
      std::cout << "spin numbering exceed nspin!" << std::endl;
    }
    if ( nspin_dim[ pos_i_lst[i] ] == 0 )
      nspin_dim[ pos_i_lst[i] ] = dim_i_lst[i];
    else if ( nspin_dim[ pos_i_lst[i] ] != dim_i_lst[i] )
    {
      std::cout << "spin dimension conflicts!" << std::endl;
    }
  }
  
  ham_dim = 1;
  for ( i = 0; i < nspin; i++ )
    ham_dim *= nspin_dim[i];
  if ( ham_dim < 2 )
  {
    std::cout << "ERROR: ham_dim = " << ham_dim << std::endl;
  }
  
  // read initial state v from data file;
  fp = fopen( vec_input_name, "rb" );
  double *vr, *vi;
  vr = new double [ ham_dim ];
  vi = new double [ ham_dim ];
  if ( fread( vr, sizeof(vr[0]), ham_dim, fp ) != ham_dim )
    std::cout << "file read error: vr" << std::endl;
  if ( fread( vi, sizeof(vi[0]), ham_dim, fp ) != ham_dim )
    std::cout << "file read error: vi" << std::endl;
  fclose( fp );
  
  v = new std::complex<double> [ ham_dim ];
  for ( i = 0; i < ham_dim; i++ )
    v[ i ] = std::complex<double>( vr[i], vi[i] );
  
  // read time line from data file;
  fp = fopen( tim_input_name, "rb" );
  size_t  tn;
  double  *ta;
  if ( fread( &tn, sizeof(tn), 1, fp ) != 1 )
    std::cout << "file read error: tn" << std::endl;
  ta = new double [ tn ];
  if ( fread( ta, sizeof(ta[0]), tn, fp ) != tn )
    std::cout << "file read error: ta" << std::endl;
  fclose( fp );
//======================================================================
  
  w = new std::complex<double> [ ham_dim ];
  
  // cuzkmv in cuda;
  cuzkmv( nspin, nTerm, coeff_lst, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, mat_i_lst, ham_dim, nspin_dim, v, pos_idx, mat_idx, w );
  
  // output to data file;
  double * w_dble;
  w_dble = new double [ ham_dim ];
  fp = fopen( output_name, "wb" );
  if ( fwrite( &ham_dim, sizeof(ham_dim), 1, fp ) != 1 )
    std::cout << "file read error: ham_dim" << std::endl;
  for ( i = 0; i < ham_dim; i++ )
    w_dble[i] = w[i].real();
  if ( fwrite( w_dble, sizeof(w_dble[0]), ham_dim, fp ) != ham_dim )
    std::cout << "file read error: w.real" << std::endl;
  for ( i = 0; i < ham_dim; i++ )
    w_dble[i] = w[i].imag();
  if ( fwrite( w_dble, sizeof(w_dble[0]), ham_dim, fp ) != ham_dim )
    std::cout << "file read error: w.imag" << std::endl;
  fclose( fp );
  
  return 0;
}

