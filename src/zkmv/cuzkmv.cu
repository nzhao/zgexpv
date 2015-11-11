/*
  Calculation of the matrix-vector multiplication exploited by CUDA, where the matrix can be decomposed into sum of the kronecker multiplication.
  
  cuzkmv( nspin, nTerm, coeff_lst, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, mat_i_lst, vlen, nspin_dim, v, w );
  
  input:
    
    nspin,                  integer, 
                            number of total spins;
    nTerm,                  integer, 
                            number of Hamiltonian terms;
    coeff_lst[ nTerm ],     double, 
                            coefficient list in each term;
    nbody_lst[ nTerm ],     size_t, 
                            spin number in each term;
    pos_idx,                size_t, 
                            total length of spin position list, the same as
                            sum of nbody_lst;
    pos_i_idx[ nTerm ],     size_t, 
                            offset of the spin position list;
    pos_i_lst[ pos_idx ],   size_t, 
                            spin position list;
    dim_i_lst[ pos_idx ],   size_t,
                            spin dimension list;
    mat_idx,                size_t, 
                            total length of the spin operator list;
    mat_i_idx[ pos_idx ],   size_t, 
                            offset of the spin operator list;
    mat_i_lst[ mat_idx ],          complex<double>, 
                            spin operator list;
    vlen,                   size_t, 
                            length of the state vector;
    v[ vlen ],              complex<double>, 
                            input state vector;
    nspin_dim[ nspin ],     size_t, 
                            list of spin dimension;
    
  output:
    
    w[ vlen ],              complex<double>, 
                            output state vector;
*/

#include <iostream>
#include <complex>
#include <cmath>

#include "cuda_runtime.h"
#include "cublas_v2.h"

void hamvec_cuda3( cublasHandle_t cublas_handle, int nspin, int nTerm, std::complex<double> *coeff_lst_zplx, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *dev_mat_i_lst, size_t vlen, size_t *nspin_dim, size_t *nspin_m_lst, size_t *nspin_n_lst, std::complex<double> *dev_v, std::complex<double> *dev_w, std::complex<double> *dev_w_med, std::complex<double> *dev_coeff_lst_zplx, size_t maxThreadsPerBlock, size_t *maxGridSize );

void cuzkmv( int nspin, int nTerm, double *coeff_lst, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *mat_i_lst, size_t vlen, size_t *nspin_dim, std::complex<double> *v, size_t pos_idx, size_t mat_idx, std::complex<double> *w )
{
  
  // coefficients for matrix;
  std::complex<double>  *coeff_lst_zplx;
  size_t                *nspin_m_lst, *nspin_n_lst;
  
  // CUDA;
  // define grid size;
  size_t maxThreadsPerBlock = 256;
  size_t maxGridSize[3]     = {2147483647,65535,65535};
  
  // variables on device;
  cudaError_t           cuda_status;
  cublasHandle_t        cublas_handle;
  std::complex<double>  *dev_mat_i_lst;
  std::complex<double>  *dev_v, *dev_w, *dev_w_med;
  std::complex<double>  *dev_coeff_lst_zplx;
  
  size_t                i, j, k;
  
  coeff_lst_zplx = new std::complex<double> [ nTerm ];
  for ( i = 0; i < nTerm; i++ )
    coeff_lst_zplx[i] = std::complex<double> ( coeff_lst[i], 0.0  );// h;
//    coeff_lst_zplx[i] = std::complex<double> ( 0.0, -coeff_lst[i] );// -i*h;
  
  nspin_m_lst = new size_t [ nspin ];
  nspin_n_lst = new size_t [ nspin ];
  for ( i = 0; i < nspin; i++ )
  {
    k = 1;
    for ( j = 0; j < i; j++ )
      k = k * nspin_dim[ j ];
    nspin_m_lst[ i ] = k;
    k = 1;
    for ( j = i + 1; j < nspin; j++ )
      k = k * nspin_dim[ j ];
    nspin_n_lst[ i ] = k;
  }
  
  // malloc on device;
  cuda_status = cudaMalloc((void**)&dev_mat_i_lst, mat_idx * sizeof( mat_i_lst[0] ) );
  if (cuda_status != cudaSuccess)
    std::cout << "Device malloc failed: dev_mat_i_lst" << std::endl;
  cuda_status = cudaMalloc((void**)&dev_v, vlen * sizeof( v[0] ) );
  if (cuda_status != cudaSuccess)
    std::cout << "Device malloc failed: dev_v" << std::endl;
  cuda_status = cudaMalloc((void**)&dev_w, vlen * sizeof( w[0] ) );
  if (cuda_status != cudaSuccess)
    std::cout << "Device malloc failed: dev_w" << std::endl;
  cuda_status = cudaMalloc((void**)&dev_w_med, vlen * sizeof( v[0] ) );
  if (cuda_status != cudaSuccess)
    std::cout << "Device malloc failed: dev_w_med" << std::endl;
  cuda_status = cudaMalloc((void**)&dev_coeff_lst_zplx, nTerm * sizeof( coeff_lst_zplx[0] ) );
  if (cuda_status != cudaSuccess)
    std::cout << "Device malloc failed: dev_coeff_lst_zplx" << std::endl;
  
  // memcpy to device;
  cuda_status = cudaMemcpy( dev_mat_i_lst, mat_i_lst, ( mat_idx * sizeof( mat_i_lst[0] ) ), cudaMemcpyHostToDevice );
  if (cuda_status != cudaSuccess)
    std::cout << "Device memcpy failed: dev_mat_i_lst" << std::endl;
  cuda_status = cudaMemcpy( dev_v, v, ( vlen * sizeof( v[0] ) ), cudaMemcpyHostToDevice );
  if (cuda_status != cudaSuccess)
    std::cout << "Device memcpy failed: dev_v" << std::endl;
  cuda_status = cudaMemcpy( dev_w, v, ( vlen * sizeof( v[0] ) ), cudaMemcpyHostToDevice );
  if (cuda_status != cudaSuccess)
    std::cout << "Device memcpy failed: dev_w" << std::endl;
  cuda_status = cudaMemcpy( dev_w_med, v, ( vlen * sizeof( v[0] ) ), cudaMemcpyHostToDevice );
  if (cuda_status != cudaSuccess)
    std::cout << "Device memcpy failed: dev_w_med" << std::endl;
  cuda_status = cudaMemcpy( dev_coeff_lst_zplx, coeff_lst_zplx, ( nTerm * sizeof( coeff_lst_zplx[0] ) ), cudaMemcpyHostToDevice );
  if (cuda_status != cudaSuccess)
    std::cout << "Device memcpy failed: dev_coeff_lst_zplx" << std::endl;
  
  cublasCreate( &cublas_handle );
  
  hamvec_cuda3( cublas_handle, nspin, nTerm, coeff_lst_zplx, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, dev_mat_i_lst, vlen, nspin_dim, nspin_m_lst, nspin_n_lst, dev_v, dev_w, dev_w_med, dev_coeff_lst_zplx, maxThreadsPerBlock, maxGridSize );
  
  cuda_status = cudaMemcpy( w, dev_w, ( vlen * sizeof( w[0]) ), cudaMemcpyDeviceToHost);
  if (cuda_status != cudaSuccess)
    std::cout << "Device memcpy failed: w" << std::endl;
  
  cublasDestroy( cublas_handle );
  cudaFree( dev_mat_i_lst );
  cudaFree( dev_v );
  cudaFree( dev_w );
  cudaFree( dev_w_med );
  cudaFree( dev_coeff_lst_zplx );
  
}
