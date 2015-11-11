      SUBROUTINE MAIN_CUDA( nspin, nTerm, coeff_lst, nbody_lst, 
     .  pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, mat_i_lst, ham_dim, 
     .  nspin_dim, v, pos_idx, mat_idx, k, maxThreadsPerBlock, 
     .  maxGridSize, tn, ta, m, tol, itrace, w_seq, w_seq_len )
      
!     Sample program to take use of cuda_zgexpv, originated from
!     sample_z in expokit. The CMPLX to DCMPLX problem should be
!     corrected in the subroutine of ZGPADM in expokit.f
*
*---  sample program illustrating the use of ZGEXPV and ZHEXPV ...
*     Hermitian problem (Example 6.2 in the Expokit report) ...
*
      implicit none
!      external zgcoov, zgcrsv, zgccsv

!      double precision tic, tac, clock

*---  matrix data ...
      INTEGER                 ::  n, nz
      
      INTEGER         w_seq_len
      COMPLEX*16  ::  w_seq(w_seq_len)

*---  arguments variables ...
      integer m, mmax, lwsp, liwsp
      parameter( mmax = 50 )
      parameter( liwsp = mmax+2 )
      integer count1, count2, count3, count4, countr
      integer iwsp(liwsp)
      double precision t, tol, anorm
      COMPLEX*16, ALLOCATABLE ::  wsp(:)
      
      ! Hamiltonian coefficients;
      INTEGER                 ::  nspin, nTerm, tn
      INTEGER*8               ::  pos_idx, mat_idx, ham_dim
      REAL*8                  ::  coeff_lst(nTerm), ta(tn)
      INTEGER*8               ::  nbody_lst(nTerm), pos_i_idx(nTerm)
      INTEGER*8               ::  pos_i_lst(pos_idx),dim_i_lst(pos_idx), 
     .                            mat_i_idx(pos_idx)
      INTEGER*8               ::  nspin_dim(nspin)
      COMPLEX*16              ::  mat_i_lst(mat_idx), v(ham_dim)
      
      ! CUDA, limit the grid and block size;
      INTEGER*8 ::  maxThreadsPerBlock
      INTEGER*8 ::  maxGridSize(3)
      
      ! Lanczos method subspace size, to get the largest eigenvalue of 
      ! the Hamiltonian;
      INTEGER                 ::  k
      
      integer i, j, itrace, iflag
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      
*
*---  Set input arguments ...
*      
      write(*,*) nspin, nTerm, ham_dim, pos_idx, mat_idx, k, tn, m,
     . w_seq_len
      n = ham_dim;
      
      lwsp = n*(mmax+2)+5*(mmax+2)**2+7
      ALLOCATE( wsp(lwsp) )
      
      if ( m.gt.mmax ) stop 'Please increase mmax.'
      
*
*---  Executable statements ...
*
      
      ! an approximation of some norm of A, similar to infinity-norm;
      ! recalculated in hamvec_zgexpv;
      anorm = 0;
      
*---  compute w = exp(t*A)v with ZGEXPV ...
      
      !=================================================================
      ! CUDA
      
      CALL SYSTEM_CLOCK( COUNT = count1, COUNT_RATE = countr )
      
      CALL hamvec_zgexpv_w_cuda_profile( nspin, nTerm, coeff_lst, 
     .  nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, 
     .  mat_i_lst, ham_dim, nspin_dim, v, pos_idx, mat_idx, k, 
     .  maxThreadsPerBlock, maxGridSize,
     .  tn, ta, m, tol, anorm,
     .  wsp, lwsp, iwsp, liwsp, itrace, iflag, w_seq, w_seq_len )
      
      CALL SYSTEM_CLOCK( COUNT = count3, COUNT_RATE = countr )
      WRITE(*,*) 'nspin = ', int(nspin,2), 
     .  ', Elapse time in zgexpv_cuda is', 
     .  REAL(DBLE(count3-count1)/DBLE(countr),4), 's'
      
      !=================================================================
      
      print 9003,'nstep     =',iwsp(4)
*---  
      print 9001,'----------------------------------------------------'
      print 9001,'ZGEXPV has completed:'
      print 9001,'----------------------------------------------------'
!      print 9001,'w(1:10) ='
!      do i = 1,10
!         print *, w(i)
!      enddo
*---  display some statistics if desired ...
      print 9001,'final report----------------------------------------'
!      print 9002,'runtime   = ',tac-tic
      print 9002,'||A||_inf = ',anorm
!      print 9003,'nz        =',nz
      print 9003,'n         =',n
      print 9003,'m         =',m
      print 9003,'itrace    =',itrace
      print 9003,'iflag     =',iflag
      print 9003,'ibrkflag  =',iwsp(6)
      print 9003,'mbrkdwn   =',iwsp(7)
      print 9003,'nstep     =',iwsp(4)
      print 9003,'nreject   =',iwsp(5)
      print 9003,'nmult     =',iwsp(1)
      print 9003,'nexph     =',iwsp(2)
      print 9003,'nscale    =',iwsp(3)
      print 9002,'tol       = ',tol
      print 9002,'t         = ',t
      print 9002,'tbrkdwn   = ',DBLE( wsp(7) )
      print 9002,'step_min  = ',DBLE( wsp(1) )
      print 9002,'step_max  = ',DBLE( wsp(2) )
      print 9002,'max_round = ',DBLE( wsp(3) )
      print 9002,'sum_round = ',DBLE( wsp(4) )
      print 9002,'max_error = ',DBLE( wsp(5) )
      print 9002,'sum_error = ',DBLE( wsp(6) )
      print 9002,'hump      = ',DBLE( wsp(9) )
      print 9002,'scale-norm= ',DBLE( wsp(10) )
      
 9001 format(A)
 9002 format(A,E9.2)
 9003 format(A,I9)
      END SUBROUTINE MAIN_CUDA
*----------------------------------------------------------------------|
