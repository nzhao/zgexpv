// |w> = |w> + H * |v>;
// H = sum(H_i) + sum(H_ij);
// H_i = h_i * S_(i,xyz) or h_ij * S_(i,xyz) * S_(j,xyz);
// w_med = H_ij * w_med 

#include <iostream>
#include <complex>
#include <cmath>

#include "cuda_runtime.h"
#include "cublas_v2.h"

// global variable used for texture memory optimization;
texture< int2, 1, cudaReadModeElementType > texRef;

//======================================================================
// kron_cuda_v1, dev_v -> dev_w_med;
//======================================================================
__global__ void kron_cuda_v1( const size_t m, const size_t s, const size_t n, const cuDoubleComplex *A, size_t mat_i_idx_idx, cuDoubleComplex *x, cuDoubleComplex *w )
{
  cuDoubleComplex res;
  cuDoubleComplex mid;
  int2 a1,a2;  

  size_t k;
  size_t q, idx2;
  
  extern __shared__ cuDoubleComplex x_shd[ ];
  
  k = blockDim.y * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.y;
  
//  idx2 = s * n * ( k / n ) + k % n;
  idx2 = ( s - 1 ) * n * ( k / n ) + k;
  
  // copy x to the shared memory x_shd;
  x_shd[ threadIdx.x + blockDim.x * threadIdx.y ] = x[ idx2 + n * threadIdx.x];
  __syncthreads();
  
  // matrix multiplication using the shared memory;
  //res = cuCmul( x_shd[ blockDim.x * threadIdx.y ], A[ threadIdx.x ] );
  a1 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x ) * 2 );
  a2 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x ) * 2 + 1 );
  mid.x = __hiloint2double( a1.y, a1.x );
  mid.y = __hiloint2double( a2.y, a2.x );
  
  res = cuCmul( x_shd[ blockDim.x * threadIdx.y ], mid);
  for ( q = 1; q < s; q++ )
  {
    a1 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x + s * q ) * 2 );
    a2 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x + s * q ) * 2 + 1 );
    mid.x = __hiloint2double( a1.y, a1.x );
    mid.y = __hiloint2double( a2.y, a2.x );
    res = cuCadd( res, cuCmul( x_shd[ q + blockDim.x * threadIdx.y ], mid) );
  }
  
  // copy to the global memory x;
  w[ idx2 + n * threadIdx.x ] = res;
}
//======================================================================

//======================================================================
// kron_cuda_v2, dev_w_med -> dev_w_med;
//======================================================================
__global__ void kron_cuda_v2( const size_t m, const size_t s, const size_t n, const cuDoubleComplex *A, size_t mat_i_idx_idx, cuDoubleComplex *x )
{
  cuDoubleComplex res;
  cuDoubleComplex mid;
  
  int2   a1, a2;  

  size_t k;
  size_t q, idx2;
  
  extern __shared__ cuDoubleComplex x_shd[ ];
//  extern __shared__ cuDoubleComplex y_shd[ ];
  
  k = blockDim.y * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.y;
  
//  idx2 = s * n * ( k / n ) + k % n;
  idx2 = ( s - 1 ) * n * ( k / n ) + k; 
  
  // copy x to the shared memory x_shd;
  x_shd[ threadIdx.x + blockDim.x * threadIdx.y ] = x[ idx2 + n * threadIdx.x];
  __syncthreads();
  
  // matrix multiplication using the shared memory;
  //res = cuCmul( x_shd[ blockDim.x * threadIdx.y ], A[ threadIdx.x ] );
  a1 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x ) * 2 );
  a2 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x ) * 2 + 1 );
  mid.x = __hiloint2double( a1.y, a1.x );
  mid.y = __hiloint2double( a2.y, a2.x );
  
  res = cuCmul( x_shd[ blockDim.x * threadIdx.y ], mid );
  for ( q = 1; q < s; q++ )
  {
    a1 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x + s * q ) * 2 );
    a2 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x + s * q ) * 2 + 1 );
    mid.x = __hiloint2double( a1.y, a1.x );
    mid.y = __hiloint2double( a2.y, a2.x );
    res = cuCadd( res, cuCmul( x_shd[ q + blockDim.x * threadIdx.y ], mid ) );
  }
//  y_shd[ threadIdx.x + blockDim.x * threadIdx.y ] = res;
//  __syncthreads();
  
  // copy y_shd to the global memory x;
//  x[ idx2 + n * threadIdx.x ] = y_shd[ threadIdx.x + blockDim.x * threadIdx.y ];
  x[ idx2 + n * threadIdx.x ] = res;
}
//======================================================================

//======================================================================
// kron_cuda_v3, dev_w_med -> dev_w;
//======================================================================
__global__ void kron_cuda_v3( const size_t m, const size_t s, const size_t n, const cuDoubleComplex *A, size_t mat_i_idx_idx, cuDoubleComplex *x, cuDoubleComplex *w, cuDoubleComplex *coeff_lst, size_t nT )
{
  cuDoubleComplex res;
  cuDoubleComplex mid;
  int2 a1,a2;  

  size_t k;
  size_t q, idx2;
  
  extern __shared__ cuDoubleComplex x_shd[ ];
  
  k = blockDim.y * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.y;
  
//  idx2 = s * n * ( k / n ) + k % n;
  idx2 = ( s - 1 ) * n * ( k / n ) + k;
  
  // copy x to the shared memory x_shd;
  x_shd[ threadIdx.x + blockDim.x * threadIdx.y ] = x[ idx2 + n * threadIdx.x];
  __syncthreads();
  
  // matrix multiplication using the shared memory;
  //res = cuCmul( x_shd[ blockDim.x * threadIdx.y ], A[ threadIdx.x ] );
  a1 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x ) * 2 );
  a2 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x ) * 2 + 1 );
  mid.x = __hiloint2double( a1.y, a1.x );
  mid.y = __hiloint2double( a2.y, a2.x );
  
  res = cuCmul( x_shd[ blockDim.x * threadIdx.y ], mid);
  for ( q = 1; q < s; q++ )
  {
    a1 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x + s * q ) * 2 );
    a2 = tex1Dfetch( texRef, ( mat_i_idx_idx + threadIdx.x + s * q ) * 2 + 1 );
    mid.x = __hiloint2double( a1.y, a1.x );
    mid.y = __hiloint2double( a2.y, a2.x );
    res = cuCadd( res, cuCmul( x_shd[ q + blockDim.x * threadIdx.y ], mid) );
  }
  
  // copy to the global memory x;
  res = cuCmul( res, coeff_lst[ nT ] );
  w[ idx2 + n * threadIdx.x ] = cuCadd( w[ idx2 + n * threadIdx.x ], res );
}
//======================================================================

__global__ void vecrzt_kernel( cuDoubleComplex *x )
{
  size_t  idx;
  
  idx = blockDim.x * ( blockIdx.x + gridDim.x * blockIdx.y + gridDim.x * gridDim.y * blockIdx.z ) + threadIdx.x;
  
  x[ idx ] = make_cuDoubleComplex(0.0,0.0);
}

void hamvec_cuda3( cublasHandle_t cublas_handle, int nspin, int nTerm, std::complex<double> *coeff_lst_zplx, size_t *nbody_lst, size_t *pos_i_idx, size_t *pos_i_lst, size_t *dim_i_lst, size_t *mat_i_idx, std::complex<double> *dev_mat_i_lst, size_t vlen, size_t *nspin_dim, size_t *nspin_m_lst, size_t *nspin_n_lst, std::complex<double> *dev_v, std::complex<double> *dev_w, std::complex<double> *dev_w_med, std::complex<double> *dev_coeff_lst_zplx, size_t maxThreadsPerBlock, size_t *maxGridSize )
{
/*
  Calculate the action of a Hamiltonian on a state vector.
  The Hamiltonian can be decomposed into many terms, where each term 
  consists of 
  
  input:
    
    cublas_handle,  handle of cublas;
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
    dev_mat_i_lst,  list of the operators of bodies in each interaction; 
    vlen,           dimension of the state vector;
    nspin_dim,      list of the dimension of each body;
    nspin_m_lst,    list of the dimension of first m bodies;
    nspin_n_lst,    list of the dimension of last n bodies;
    dev_v,          input state vector;
    dev_w_med,      intermediate state vector;
    dev_coeff_lst_zplx, the same as "coeff_lst_zplx";
    maxThreadsPerBlock,   max block size for CUDA;
    maxGridSize,    grid size limit for CUDA;
  
  output:
    
    dev_w,          output state vector, after the Hamiltonian acting on
                    the input state;
*/
  
  size_t  nT, nbody, nb;
  size_t  idx, pos_i, dim_i;
  std::complex<double>  coeff;
  
  size_t  m, n;
  
  dim3    grid_dim, block_dim;
  size_t  dimex1, dimex2;
  
  // optimization by texture memory, bind the operator list;
  cudaBindTexture( 0, texRef, dev_mat_i_lst );
  
  //====================================================================
  // initialization of w;
  //====================================================================
  
  // set block size, 2D;
  block_dim.z = 1;
  block_dim.y = 1;
  dimex1 = maxThreadsPerBlock;
  dimex2 = vlen;
  while ( (dimex2 % dimex1) != 0 ) dimex1--;
  block_dim.x = dimex1;
  dimex2 /= dimex1;
  if ( dimex2 <= maxGridSize[0] * maxGridSize[1] * maxGridSize[2] )
  {
    dimex1 = dimex2 / ( maxGridSize[0] * maxGridSize[1] );
    if ( ( dimex2 % ( maxGridSize[0] * maxGridSize[1] ) ) > 0 )
      dimex1++;
    while ( (dimex2 % dimex1) != 0 ) dimex1++;
    grid_dim.z = dimex1;
    dimex2 = dimex2 / grid_dim.z;
    dimex1 = dimex2 / maxGridSize[0];
    if ( ( dimex2 % maxGridSize[0] ) > 0)
      dimex1++;
    while ( (dimex2 % dimex1) != 0 ) dimex1++;
    grid_dim.y = dimex1;
    grid_dim.x = dimex2 / grid_dim.y;
  }
  else
  {
    std::cout << "block number exceeds limit." << std::endl;
    return;
  }
  
  // reset the output vector;
  vecrzt_kernel<<< grid_dim, block_dim >>>( (cuDoubleComplex*)dev_w );
  
  //====================================================================
  // action of the Hamiltonian on the input state;
  //====================================================================
  
  // loop over each term of the Hamiltonian;
  for ( nT = 0; nT < nTerm; nT++ )
  {
    coeff = coeff_lst_zplx[ nT ];
    nbody = nbody_lst[ nT ];
    
    if ( abs(coeff) == 0 ) continue;
    
    // loop over each body in each term of the Hamiltonian;
    for ( nb = 0; nb < nbody; nb++ )
    {
      idx = pos_i_idx[ nT ] + nbody - 1 - nb;
      pos_i = pos_i_lst[ idx ];
      dim_i = dim_i_lst[ idx ];
      
      m = nspin_m_lst[ pos_i ];
      n = nspin_n_lst[ pos_i ];
      
      // set block size, 2D;
      block_dim.z = 1;
      block_dim.x = dim_i;
      dimex1 = maxThreadsPerBlock / block_dim.x;
      dimex2 = m * n;
      while ( (dimex2 % dimex1) != 0 ) dimex1--;
      block_dim.y = dimex1;
      // set grid size, 1D-3D;
      dimex2 /= dimex1;// grid size is m * n / block_dim.y;
      if ( dimex2 <= maxGridSize[0] * maxGridSize[1] * maxGridSize[2] )
      {
        dimex1 = dimex2 / ( maxGridSize[0] * maxGridSize[1] );
        if ( ( dimex2 % ( maxGridSize[0] * maxGridSize[1] ) ) > 0 )
          dimex1++;
        while ( (dimex2 % dimex1) != 0 ) dimex1++;
        grid_dim.z = dimex1;
        dimex2 = dimex2 / grid_dim.z;
        dimex1 = dimex2 / maxGridSize[0];
        if ( ( dimex2 % maxGridSize[0] ) > 0)
          dimex1++;
        while ( (dimex2 % dimex1) != 0 ) dimex1++;
        grid_dim.y = dimex1;
        grid_dim.x = dimex2 / grid_dim.y;
      }
      else
      {
        std::cout << "block number exceeds limit." << std::endl;
        return;
      }
      
      if ( nbody == 1 )
        kron_cuda_v3<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_v, (cuDoubleComplex*)dev_w, (cuDoubleComplex*)dev_coeff_lst_zplx, nT );
      else
      {
        if ( nb == 0 )
          kron_cuda_v1<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_v, (cuDoubleComplex*)dev_w_med );
        else if ( nb == nbody - 1 )
          kron_cuda_v3<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_w_med, (cuDoubleComplex*)dev_w, (cuDoubleComplex*)dev_coeff_lst_zplx, nT );
        else
          kron_cuda_v2<<< grid_dim, block_dim, block_dim.x * block_dim.y * sizeof(std::complex<double>) >>>( m, dim_i, n, (cuDoubleComplex*)dev_mat_i_lst, mat_i_idx[idx], (cuDoubleComplex*)dev_w_med );
      }
      
    }
    
  }
  
  // optimization by texture memory, release the bind;
  cudaUnbindTexture( texRef );
}

//======================================================================
// interface for FORTRAN;
//======================================================================

// declaration of function handle;
#define HAMVEC_CUDA3       hamvec_cuda3_
#define HAMVEC_CUDA3_INIT  hamvec_cuda3_init_
#define HAMVEC_CUDA3_TERM  hamvec_cuda3_term_

// delcaration of function;
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */

void HAMVEC_CUDA3( size_t *cublas_handle_ptr, int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, size_t *dev_mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, size_t *dev_v_ptr, size_t *dev_w_ptr, size_t *dev_w_med_ptr, size_t *dev_coeff_lst_zplx_ptr, size_t *maxThreadsPerBlock_ptr, size_t *maxGridSize_ptr );

void HAMVEC_CUDA3_INIT( size_t *cublas_handle_ptr );

void HAMVEC_CUDA3_TERM( size_t *cublas_handle_ptr );

#if defined(__cplusplus)
}
#endif /* __cplusplus */

// interface of function;
void HAMVEC_CUDA3( size_t *cublas_handle_ptr, int *nspin_ptr, int *nTerm_ptr, std::complex<double> *coeff_lst_zplx_ptr, size_t *nbody_lst_ptr, size_t *pos_i_idx_ptr, size_t *pos_i_lst_ptr, size_t *dim_i_lst_ptr, size_t *mat_i_idx_ptr, size_t *dev_mat_i_lst_ptr, size_t *ham_dim_ptr, size_t *nspin_dim_ptr, size_t *nspin_m_lst_ptr, size_t *nspin_n_lst_ptr, size_t *dev_v_ptr, size_t *dev_w_ptr, size_t *dev_w_med_ptr, size_t *dev_coeff_lst_zplx_ptr, size_t *maxThreadsPerBlock_ptr, size_t *maxGridSize_ptr )
{
  cublasHandle_t        cublas_handle   = (cublasHandle_t)*cublas_handle_ptr;
  int                   nspin           = *nspin_ptr;
  int                   nTerm           = *nTerm_ptr;
  std::complex<double>  *coeff_lst_zplx = coeff_lst_zplx_ptr;
  size_t                *nbody_lst      = nbody_lst_ptr;
  size_t                *pos_i_idx      = pos_i_idx_ptr;
  size_t                *pos_i_lst      = pos_i_lst_ptr;
  size_t                *dim_i_lst      = dim_i_lst_ptr;
  size_t                *mat_i_idx      = mat_i_idx_ptr;
  std::complex<double>  *dev_mat_i_lst  = (std::complex<double>*)(*dev_mat_i_lst_ptr);
  size_t                vlen            = *ham_dim_ptr;
  size_t                *nspin_dim      = nspin_dim_ptr;
  size_t                *nspin_m_lst    = nspin_m_lst_ptr;
  size_t                *nspin_n_lst    = nspin_n_lst_ptr;
  std::complex<double>  *dev_v          = (std::complex<double>*)(*dev_v_ptr);
  std::complex<double>  *dev_w          = (std::complex<double>*)(*dev_w_ptr);
  std::complex<double>  *dev_w_med      = (std::complex<double>*)(*dev_w_med_ptr);
  std::complex<double>  *dev_coeff_lst_zplx = (std::complex<double>*)(*dev_coeff_lst_zplx_ptr);
  size_t            maxThreadsPerBlock  = *maxThreadsPerBlock_ptr;
  size_t                *maxGridSize    = maxGridSize_ptr;
  
  hamvec_cuda3( cublas_handle, nspin, nTerm, coeff_lst_zplx, nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, dev_mat_i_lst, vlen, nspin_dim, nspin_m_lst, nspin_n_lst, dev_v, dev_w, dev_w_med, dev_coeff_lst_zplx, maxThreadsPerBlock, maxGridSize );
}

void HAMVEC_CUDA3_INIT( size_t *cublas_handle_ptr )
{
  // initialization of cublas handle for FORTRAN;
  cublasHandle_t cublas_handle;
  cublasCreate( &cublas_handle );
  *cublas_handle_ptr = (size_t)cublas_handle;
}

void HAMVEC_CUDA3_TERM( size_t *cublas_handle_ptr )
{
  // termination of cublas handle for FORTRAN;
  cublasHandle_t cublas_handle;
  cublas_handle = (cublasHandle_t)*cublas_handle_ptr;
  cublasDestroy( cublas_handle );
}

