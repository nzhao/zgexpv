// |w> = |w> + H * |v>;
// H = sum(H_i) + sum(H_ij);
// H_i = h_i * S_(i,xyz) or h_ij * S_(i,xyz) * S_(j,xyz);
// w_med = H_ij * w_med 

#include <iostream>
#include <complex>
#include <cmath>

#include "omp.h"
#define MKL_Complex16 std::complex<double>
#include "mkl.h"

//======================================================================
// kron_func_v2.cpp
//======================================================================
void kron_func_v2( size_t m, size_t s, size_t n, std::complex<double> *A, std::complex<double> *x )
{
  std::complex<double> zero(0.0,0.0);
  std::complex<double> one(1.0,0.0);
  
  size_t inc1, inc2, inc3;
  
  inc1 = s * n;
  inc2 = n;
  inc3 = s;
  
  #pragma omp parallel
  {
  size_t k;
  size_t i, j, p, q, idx1, idx2;
  std::complex<double> res;
  std::complex<double> y[s];
  #pragma omp for
  for ( k = 0; k < m*n; k++ )
  {
    i = k / n;
    j = k % n;
//  for ( i = 0; i < m; i++ )
//  {
    idx1 = inc1 * i;
//    for ( j = 0; j < n; j++ )
//    {
      idx2 = idx1 + j;
      
      for ( p = 0; p < s; p++ )
      {
//        idx3 = s * p;
        res = zero;
        for ( q = 0; q < s; q++ )
        {
          res += x[ idx2 + inc2 * q ] * A[ p + inc3 * q ];
        }
        y[ p ] = res;
      }
      for ( p = 0; p < s; p++ )
        x[ idx2 + inc2 * p] = y[p];
      
//    }
//  }
  }
  } //  end of pragma omp parallel;
  
}
//======================================================================

void hamvec_func3( int nspin, int nTerm, std::complex<double> *coeff_lst_zplx, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *mat_i_lst, size_t vlen, size_t *nspin_dim, size_t *nspin_m_lst, size_t *nspin_n_lst, std::complex<double> *v, std::complex<double> *w, std::complex<double> *w_med )
{
/*
  Calculate the action of a Hamiltonian on a state vector.
  The Hamiltonian can be decomposed into many terms, where each term 
  consists of 
  
  input:
    
    nspin,          number of bodies;
    nTerm,          number of interactions in the Hamiltonian;
    coeff_lst_zplx, list of the coefficient in each interaction;
    nbody_lst,      list of the number of bodies in each interaction;
    pos_i_idx,      list of the index of the position list of the body 
                    in each interaction;
    pos_i_lst,      list of the positions of bodies in each interaction;
    dim_i_lst,      list of the dimension of operator of each body in 
                    each interaction;
    mat_i_idx,      list of the index of the operator list of the body
                    in each interaction;
    mat_i_lst,      list of the operators of bodies in each interaction; 
    vlen,           dimension of the state vector;
    nspin_dim,      list of the dimension of each body;
    nspin_m_lst,    list of the dimension of first m bodies;
    nspin_n_lst,    list of the dimension of last n bodies;
    v,              input state vector;
    w_med,          intermediate state vector;
  
  output:
    
    w,              output state vector, after the Hamiltonian acting on
                    the input state;
*/
  size_t                nT, nbody, nb;
  size_t                idx, pos_i, dim_i;
  std::complex<double>  coeff;
  
  size_t                i, m, n;
  
  // initialization of w;
  #pragma omp parallel
  {
  #pragma omp for
  for ( i = 0; i < vlen; i++ )
    w[i] = std::complex<double>(0.0,0.0);
  }
  
  for ( nT = 0; nT < nTerm; nT++ )
  {
    coeff = coeff_lst_zplx[ nT ];
    nbody = nbody_lst[ nT ];
    
    cblas_zcopy( vlen, v, 1, w_med, 1 );
    for ( nb = 0; nb < nbody; nb++ )
    {
      idx = pos_i_idx[ nT ] + nbody - 1 - nb;
      pos_i = pos_i_lst[ idx ];
      dim_i = dim_i_lst[ idx ];
      
      m = nspin_m_lst[ pos_i ];
      n = nspin_n_lst[ pos_i ];
      
      kron_func_v2( m, dim_i, n, &mat_i_lst[ mat_i_idx[ idx ] ], w_med );
    }
    cblas_zaxpy( vlen, &coeff, w_med, 1, w, 1 );
  }
}

// interface for FORTRAN;

#define HAMVEC_FUNC3       hamvec_func3_

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

void HAMVEC_FUNC3( int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, std::complex<double> *mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, std::complex<double> *v_ptr, std::complex<double> *w_ptr, std::complex<double> *w_med_ptr );

#if defined(__cplusplus)
}
#endif /* __cplusplus */

void HAMVEC_FUNC3( int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, std::complex<double> *mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, std::complex<double> *v_ptr, std::complex<double> *w_ptr, std::complex<double> *w_med_ptr )
{
  int                   nspin           = *nspin_ptr;
  int                   nTerm           = *nTerm_ptr;
  std::complex<double>  *coeff_lst_zplx = coeff_lst_zplx_ptr;
  size_t                *nbody_lst      = nbody_lst_ptr;
  size_t                *pos_i_idx      = pos_i_idx_ptr;
  size_t                *pos_i_lst      = pos_i_lst_ptr;
  size_t                *dim_i_lst      = dim_i_lst_ptr;
  size_t                *mat_i_idx      = mat_i_idx_ptr;
  std::complex<double>  *mat_i_lst      = mat_i_lst_ptr;
  size_t                vlen            = *ham_dim_ptr;
  size_t                *nspin_dim      = nspin_dim_ptr;
  size_t                *nspin_m_lst    = nspin_m_lst_ptr;
  size_t                *nspin_n_lst    = nspin_n_lst_ptr;
  std::complex<double>  *v              = v_ptr;
  std::complex<double>  *w              = w_ptr;
  std::complex<double>  *w_med          = w_med_ptr;
  
  hamvec_func3( nspin, nTerm, coeff_lst_zplx, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, mat_i_lst, vlen, nspin_dim, nspin_m_lst, nspin_n_lst, v, w, w_med );
}

