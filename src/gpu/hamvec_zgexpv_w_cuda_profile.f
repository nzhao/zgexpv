*----------------------------------------------------------------------|
*	 subroutine ZGEXPV( n, m, t, v, w, tol, anorm,
*     .                   wsp,lwsp, iwsp,liwsp, matvec, itrace,iflag )
      SUBROUTINE HAMVEC_ZGEXPV_W_CUDA_PROFILE( nspin, nTerm, coeff_lst, 
     .  nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, 
     .  mat_i_lst, ham_dim, nspin_dim, v, pos_idx, mat_idx, k, 
     .  maxThreadsPerBlock, maxGridSize,
     .  tn, ta, m, tol,
     .  anorm, wsp, lwsp, iwsp, liwsp, itrace, iflag, w_seq, w_seq_len )

      implicit none
      integer          n, nz
      INTEGER           tn
      REAL*8            ta(tn)
	  
      integer          m, lwsp, liwsp, itrace, iflag, iwsp(liwsp)
      double precision t, tol, anorm
      complex*16       wsp(lwsp)
      INTEGER         w_seq_len
      COMPLEX*16      w_seq(w_seq_len)
      INTEGER           tt
      
*      external         matvec
*-----Purpose----------------------------------------------------------|
*
*---  ZGEXPV computes w = exp(t*A)*v
*     for a Zomplex (i.e., complex double precision) matrix A 
*
*     It does not compute the matrix exponential in isolation but
*     instead, it computes directly the action of the exponential
*     operator on the operand vector. This way of doing so allows 
*     for addressing large sparse problems. 
*
*     The method used is based on Krylov subspace projection
*     techniques and the matrix under consideration interacts only
*     via the external routine `matvec' performing the matrix-vector 
*     product (matrix-free method).
*
*-----Arguments--------------------------------------------------------|
*
*     n      : (input) order of the principal matrix A.
*
*     m      : (input) maximum size for the Krylov basis.
*
*     t      : (input) time at wich the solution is needed (can be < 0).
*
*     v(n)   : (input) given operand vector.
*
*     w(n)   : (output) computed approximation of exp(t*A)*v.
*
*     tol    : (input/output) the requested accuracy tolerance on w. 
*              If on input tol=0.0d0 or tol is too small (tol.le.eps)
*              the internal value sqrt(eps) is used, and tol is set to
*              sqrt(eps) on output (`eps' denotes the machine epsilon).
*              (`Happy breakdown' is assumed if h(j+1,j) .le. anorm*tol)
*
*     anorm  : (input) an approximation of some norm of A.
*
*   wsp(lwsp): (workspace) lwsp .ge. n*(m+1)+n+(m+2)^2+4*(m+2)^2+ideg+1
*                                   +---------+-------+---------------+
*              (actually, ideg=6)        V        H      wsp for PADE
*                   
* iwsp(liwsp): (workspace) liwsp .ge. m+2
*
*     matvec : external subroutine for matrix-vector multiplication.
*              synopsis: matvec( x, y )
*                        complex*16 x(*), y(*)
*              computes: y(1:n) <- A*x(1:n)
*                        where A is the principal matrix.
*
*     itrace : (input) running mode. 0=silent, 1=print step-by-step info
*
*     iflag  : (output) exit flag.
*              <0 - bad input arguments 
*               0 - no problem
*               1 - maximum number of steps reached without convergence
*               2 - requested tolerance was too high
*
*-----Accounts on the computation--------------------------------------|
*     Upon exit, an interested user may retrieve accounts on the 
*     computations. They are located in the workspace arrays wsp and 
*     iwsp as indicated below: 
*
*     location  mnemonic                 description
*     -----------------------------------------------------------------|
*     iwsp(1) = nmult, number of matrix-vector multiplications used
*     iwsp(2) = nexph, number of Hessenberg matrix exponential evaluated
*     iwsp(3) = nscale, number of repeated squaring involved in Pade
*     iwsp(4) = nstep, number of integration steps used up to completion 
*     iwsp(5) = nreject, number of rejected step-sizes
*     iwsp(6) = ibrkflag, set to 1 if `happy breakdown' and 0 otherwise
*     iwsp(7) = mbrkdwn, if `happy brkdown', basis-size when it occured
*     -----------------------------------------------------------------|
*     wsp(1)  = step_min, minimum step-size used during integration
*     wsp(2)  = step_max, maximum step-size used during integration
*     wsp(3)  = x_round, maximum among all roundoff errors (lower bound) 
*     wsp(4)  = s_round, sum of roundoff errors (lower bound)
*     wsp(5)  = x_error, maximum among all local truncation errors
*     wsp(6)  = s_error, global sum of local truncation errors
*     wsp(7)  = tbrkdwn, if `happy breakdown', time when it occured
*     wsp(8)  = t_now, integration domain successfully covered
*     wsp(9)  = hump, i.e., max||exp(sA)||, s in [0,t] (or [t,0] if t<0)
*     wsp(10) = ||w||/||v||, scaled norm of the solution w.
*     -----------------------------------------------------------------|
*     The `hump' is a measure of the conditioning of the problem. The
*     matrix exponential is well-conditioned if hump = 1, whereas it is
*     poorly-conditioned if hump >> 1. However the solution can still be
*     relatively fairly accurate even when the hump is large (the hump 
*     is an upper bound), especially when the hump and the scaled norm
*     of w [this is also computed and returned in wsp(10)] are of the 
*     same order of magnitude (further details in reference below).
*
*----------------------------------------------------------------------|
*-----The following parameters may also be adjusted herein-------------|
*
      integer mxstep, mxreject, ideg
      double precision delta, gamma
      parameter( mxstep   = 0,
     .           mxreject = 0,
     .           ideg     = 6,
     .           delta    = 1.2d0,
     .           gamma    = 0.9d0 )

*     mxstep  : maximum allowable number of integration steps.
*               The value 0 means an infinite number of steps.
* 
*     mxreject: maximum allowable number of rejections at each step. 
*               The value 0 means an infinite number of rejections.
*
*     ideg    : the Pade approximation of type (ideg,ideg) is used as 
*               an approximation to exp(H). The value 0 switches to the
*               uniform rational Chebyshev approximation of type (14,14)
*
*     delta   : local truncation error `safety factor'
*
*     gamma   : stepsize `shrinking factor'
*
*----------------------------------------------------------------------|
*     Roger B. Sidje (rbs@maths.uq.edu.au)
*     EXPOKIT: Software Package for Computing Matrix Exponentials.
*     ACM - Transactions On Mathematical Software, 24(1):130-156, 1998
*----------------------------------------------------------------------|

      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )
      integer i, j, k1, mh, mx, iv, ih, j1v, ns, ifree, lfree, iexph,
     .        ireject,ibrkflag,mbrkdwn, nmult, nreject, nexph, nscale,
     .        nstep
      double precision sgn, t_out, tbrkdwn, step_min,step_max, err_loc,
     .                 s_error, x_error, t_now, t_new, t_step, t_old,
     .                 xm, beta, break_tol, p1, p2, p3, eps, rndoff,
     .                 vnorm, avnorm, hj1j, hump, SQR1
      complex*16 hij

      intrinsic AINT,ABS,CMPLX,DBLE,INT,LOG10,MAX,MIN,NINT,SIGN,SQRT
      
      ! cuda, malloc size, byte length of zomplex and integer;
      INTEGER, PARAMETER  ::  szpx = 16, sint = 4
      ! cuda funciton handle
      INTEGER     cuda_malloc
      INTEGER     cuda_memcpy_fort2c_real
      INTEGER     cuda_memcpy_fort2c_int
      INTEGER     cusparse_xcoo2csr
      INTEGER     cusparse_zcsrmv
      INTEGER     cuda_memcpy_c2fort_real
      REAL*8      cublas_dznrm2
      
      ! cuda variables
      INTEGER     cuda_status
      INTEGER*8   cusparse_handle, cusparse_descrA
      INTEGER     cusparse_create
      INTEGER     cusparse_create_mat_descr
      INTEGER     cusparse_set_mat_type
      INTEGER     cusparse_set_mat_index_base
      EXTERNAL    cusparse_destroy
      EXTERNAL    cusparse_destroy_mat_descr
      EXTERNAL    cublas_init
      EXTERNAL    cublas_shutdown
      INTEGER*8   dev_w, dev_v(m+2)
      
      ! hamvec_cuda3;
      INTEGER     nspin, nTerm
      INTEGER*8   pos_idx, mat_idx
      REAL*8      coeff_lst(nTerm)
      INTEGER*8   nbody_lst(nTerm)
      INTEGER*8   pos_i_idx(nTerm)
      INTEGER*8   pos_i_lst(pos_idx)
      INTEGER*8   dim_i_lst(pos_idx)
      INTEGER*8   mat_i_idx(pos_idx)
      COMPLEX*16  mat_i_lst(mat_idx)
      INTEGER*8   ham_dim
      INTEGER*8   nspin_dim(nspin)
      COMPLEX*16  v(ham_dim)
      COMPLEX*16, ALLOCATABLE ::  w(:), w_med(:)
      INTEGER*8   cublas_handle
      INTEGER*8   dev_w_med, dev_mat_i_lst, dev_coeff_lst_zplx
      INTEGER*8   maxThreadsPerBlock, maxGridSize(3)
      INTEGER*8,  ALLOCATABLE ::  nspin_m_lst(:), nspin_n_lst(:)
      
      COMPLEX*16, ALLOCATABLE ::  coeff_lst_zplx(:)
      
      ! Lanczos algorithm;
      INTEGER     k
      COMPLEX*16  aj, ajm
      REAL*8      bj, bjm
      REAL*8      cj, cjm
      COMPLEX*16, ALLOCATABLE ::  qj(:), qjm(:)
      REAL*8      DLARAN, DDOT, DNRM2
      INTEGER                 ::  iseed(4) = (/ 0, 0, 0, 5 /)
      REAL*8,     ALLOCATABLE ::  tk(:), qjr(:), wjr(:)
      INTEGER*8   i8
      
      REAL*8      cublas_ddot, cublas_dnrm2
      
      n = ham_dim
      nz = 10
      ALLOCATE( w(ham_dim), w_med(ham_dim) )
      
      ALLOCATE( coeff_lst_zplx(nTerm) )
      DO i = 1, nTerm
        coeff_lst_zplx(i) = DCMPLX( +coeff_lst(i), 0.0D0 )! exp(h*t)
!        coeff_lst_zplx(i) = DCMPLX( 0.0D0, -coeff_lst(i) )! exp(-i*h*t)
      END DO
      
*
*---  initialization of cuda ...
*
      ! cuda initialize function handle
      CALL hamvec_cuda3_init( cublas_handle )
      
      ! cuda, allocate memory and copy data from host to device;
      cuda_status = cuda_malloc(dev_w, n * szpx)
      DO i = 1, m + 2
        cuda_status = cuda_malloc(dev_v(i), n * szpx)
      END DO
      
      ! hamvec_cuda3;
      cuda_status = cuda_malloc(dev_w_med, ham_dim * szpx)
      cuda_status = cuda_malloc(dev_mat_i_lst, mat_idx * szpx)
      cuda_status = cuda_malloc(dev_coeff_lst_zplx, nTerm * szpx)
      
      ALLOCATE( nspin_m_lst(nspin), nspin_n_lst(nspin) )
      DO i = 1, nspin
        i8 = 1
        DO j = 1, i-1
          i8 = i8 * nspin_dim( j )
        END DO
        nspin_m_lst( i ) = i8
        i8 = 1
        DO j = i+1, nspin
          i8 = i8 * nspin_dim( j )
        END DO
        nspin_n_lst( i ) = i8
      END DO
      
      cuda_status = cuda_memcpy_fort2c_real(dev_w, v, ham_dim * szpx, 1)
      cuda_status = cuda_memcpy_fort2c_real(dev_mat_i_lst, mat_i_lst, 
     .  mat_idx * szpx, 1)
      cuda_status = cuda_memcpy_fort2c_real( dev_coeff_lst_zplx, 
     .  coeff_lst_zplx, nTerm * szpx, 1)
      
      cuda_status = cuda_memcpy_fort2c_real(dev_w_med, v, 
     .  ham_dim * szpx, 1)
      DO i = 1, m+2
        cuda_status = cuda_memcpy_fort2c_real(dev_v(i), v, 
     .  ham_dim * szpx, 1)
      END DO
      
      ! Lanczos factorization of length k; without reorthogonalization;
      IF ( k > 0 ) k = MIN( k, ham_dim )
      
      ALLOCATE( qj(ham_dim), qjm(ham_dim) )
      ALLOCATE( tk(k**2) )
      DO j = 1, k**2
        tk(j) = 0.0D0
      END DO
      DO i8 = 1, ham_dim
        qj(i8) = DCMPLX( DLARAN(iseed), DLARAN(iseed) )
        qjm(i8) = DCMPLX(0.0D0,0.0D0)
      END DO
      cuda_status = cuda_memcpy_fort2c_real( dev_v(1), qj, 
     .  ham_dim * szpx, 1)
      bj = cublas_dznrm2( ham_dim, dev_v(1), 1 )
      CALL cublas_zdscal( ham_dim, 1.0D0/bj, dev_v(1), 1 )
      cuda_status = cuda_memcpy_fort2c_real( dev_v(2), qjm, 
     .  ham_dim * szpx, 1)
      bjm = 0.0D0
      DO j = 1, k
        CALL hamvec_cuda3( cublas_handle, nspin, nTerm, coeff_lst_zplx, 
     .    nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, 
     .    dev_mat_i_lst, ham_dim, nspin_dim, nspin_m_lst, nspin_n_lst, 
     .    dev_v(1), dev_v(3), dev_w_med, dev_coeff_lst_zplx,
     .    maxThreadsPerBlock, maxGridSize )
        CALL cublas_zdotc( aj, ham_dim, dev_v(1), 1, dev_v(3), 1 )
        tk( (j-1)*k+j ) = DREAL(aj);
        CALL cublas_zaxpy( ham_dim, -aj, dev_v(1), 1, dev_v(3), 1 )
        CALL cublas_zaxpy( ham_dim, -bjm, dev_v(2), 1, dev_v(3), 1 )
        bj = cublas_dznrm2( ham_dim, dev_v(3), 1 )
        
        IF ( j .EQ. k ) EXIT
        tk( (j-1)*k+j+1) = bj
        tk( j*k+j ) = bj
        IF ( bj .LT. tol ) THEN
          EXIT
        ELSE
          CALL cublas_zcopy( ham_dim, dev_v(1), 1, dev_v(2), 1 )
          CALL cublas_zcopy( ham_dim, dev_v(3), 1, dev_v(1), 1 )
          CALL cublas_zdscal( ham_dim, 1.0D0/bj, dev_v(1), 1 )
          bjm = bj
        END IF
      END DO
      ! Power iteration;
      ALLOCATE( qjr(k), wjr(k) )
      DO j = 1, k
        qjr(j) = DLARAN(iseed)
        wjr(j) = 0.0D0
      END DO
      cjm = 0.0D0
      DO
        CALL DGEMV( 'N', k, k, 1.0D0, tk, k, qjr, 1, 0.0D0, wjr, 1 )
        cj = DDOT( k, qjr, 1, wjr, 1 )
        IF ( ABS( 1.0D0 - cjm / cj ) .LT. tol ) EXIT
        cjm = cj
        cj = DNRM2( k, wjr, 1 )
        CALL DCOPY( k, wjr, 1, qjr, 1 )
        CALL DSCAL( k, 1.0D0/cj, qjr, 1 )
      END DO
      ! cj, qjr are the eigenvalue and eigenvector;
      ! give the upper bound of the norm of the matrix;
      anorm = ABS(cj) + ABS( bj * qjr(k) )
      
      DO i = 1, nTerm
!        coeff_lst_zplx(i) = DCMPLX( +coeff_lst(i), 0.0D0 )! exp(h*t)
        coeff_lst_zplx(i) = DCMPLX( 0.0D0, -coeff_lst(i) )! exp(-i*h*t)
      END DO
      cuda_status = cuda_memcpy_fort2c_real( dev_coeff_lst_zplx, 
     .  coeff_lst_zplx, nTerm * szpx, 1)
      
*
*---  check restrictions on input parameters ...
*
      iflag = 0
      if ( lwsp.lt.n*(m+2)+5*(m+2)**2+ideg+1 ) iflag = -1
      if ( liwsp.lt.m+2 ) iflag = -2
      if ( m.ge.n .or. m.le.0 ) iflag = -3
      if ( iflag.ne.0 ) stop 'bad sizes (in input of ZGEXPV)'
*
*---  initialisations ...
*
      k1 = 2
      mh = m + 2
      iv = 1
      ih = iv + n*(m+1) + n
      ifree = ih + mh*mh
      lfree = lwsp - ifree + 1

      ibrkflag = 0
      mbrkdwn  = m
      nmult    = 0
      nreject  = 0
      nexph    = 0
      nscale   = 0

      t_out    = ABS( t )
      tbrkdwn  = 0.0d0
      step_min = t_out
      step_max = 0.0d0
      nstep    = 0
      s_error  = 0.0d0
      x_error  = 0.0d0
      t_now    = 0.0d0
      t_new    = 0.0d0

      p1 = 4.0d0/3.0d0
 1    p2 = p1 - 1.0d0
      p3 = p2 + p2 + p2
      eps = ABS( p3-1.0d0 )
      if ( eps.eq.0.0d0 ) go to 1
      if ( tol.le.eps ) tol = SQRT( eps )
      rndoff = eps*anorm

      break_tol = 1.0d-7
*>>>  break_tol = tol
*>>>  break_tol = anorm*tol
      
      sgn = SIGN( 1.0d0,t )
      
      cuda_status = cuda_memcpy_fort2c_real(dev_w, v, n * szpx, 1)
      beta = cublas_dznrm2(n, dev_w, 1)
      vnorm = beta
      hump = beta 
*
*---  obtain the very first stepsize ...
*
      SQR1 = SQRT( 0.1d0 )
      xm = 1.0d0/DBLE( m )
      p2 = tol*(((m+1)/2.72D0)**(m+1))*SQRT(2.0D0*3.14D0*(m+1))
      t_new = (1.0d0/anorm)*(p2/(4.0d0*beta*anorm))**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1
*
*---  step-by-step integration ...
*
      ! output w(t) for t in ta
      ta_output:  DO tt = 1, tn
      t_out = ta(tt)
      
 100  if ( t_now.ge.t_out ) goto 500

      nstep = nstep + 1
      t_step = MIN( t_out-t_now, t_new )
      p1 = 1.0d0/beta
      CALL cublas_zcopy(n, dev_w, 1, dev_v(1), 1)
      CALL cublas_zdscal(n, p1, dev_v(1), 1)
      do i = 1,mh*mh
         wsp(ih+i-1) = ZERO
      enddo
*
*---  Arnoldi loop ...
*
      j1v = iv + n
      do 200 j = 1,m
         nmult = nmult + 1
*         call matvec( wsp(j1v-n), wsp(j1v) )
        CALL hamvec_cuda3( cublas_handle, nspin, nTerm, coeff_lst_zplx, 
     .    nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, 
     .    dev_mat_i_lst, ham_dim, nspin_dim, nspin_m_lst, nspin_n_lst, 
     .    dev_v(j), dev_v(j+1), dev_w_med, dev_coeff_lst_zplx, 
     .    maxThreadsPerBlock, maxGridSize )
         do i = 1,j
            CALL cublas_zdotc( hij, ham_dim, dev_v(i), 1, dev_v(j+1), 1)
            CALL cublas_zaxpy(ham_dim, -hij, dev_v(i), 1, dev_v(j+1), 1)
            wsp(ih+(j-1)*mh+i-1) = hij
         enddo
         hj1j = cublas_dznrm2( ham_dim, dev_v(j+1), 1 )
      
*---     if `happy breakdown' go straightforward at the end ... 
         if ( hj1j.le.break_tol ) then
            print*,'happy breakdown: mbrkdwn =',j,' h =',hj1j
            k1 = 0
            ibrkflag = 1
            mbrkdwn = j
            tbrkdwn = t_now
            t_step = t_out-t_now
            goto 300
         endif
         wsp(ih+(j-1)*mh+j) = DCMPLX( hj1j )
         
         CALL cublas_zdscal(n, 1.0d0/hj1j, dev_v(j+1), 1)

         j1v = j1v + n
 200  continue
      nmult = nmult + 1
*      call matvec( wsp(j1v-n), wsp(j1v) )
      CALL hamvec_cuda3( cublas_handle, nspin, nTerm, coeff_lst_zplx, 
     .    nbody_lst, pos_i_idx, pos_i_lst, dim_i_lst, mat_i_idx, 
     .    dev_mat_i_lst, ham_dim, nspin_dim, nspin_m_lst, nspin_n_lst, 
     .    dev_v(m+1), dev_v(m+2), dev_w_med, dev_coeff_lst_zplx, 
     .    maxThreadsPerBlock, maxGridSize )
      avnorm = cublas_dznrm2( ham_dim, dev_v(m+2), 1 )
      
*
*---  set 1 for the 2-corrected scheme ...
*
 300  continue
      wsp(ih+m*mh+m+1) = ONE
*
*---  loop while ireject<mxreject until the tolerance is reached ...
*
      ireject = 0
 401  continue
*
*---  compute w = beta*V*exp(t_step*H)*e1 ...
*
      nexph = nexph + 1
      mx = mbrkdwn + k1
      if ( ideg.ne.0 ) then
*---     irreducible rational Pade approximation ...
         call ZGPADM( ideg, mx, sgn*t_step, wsp(ih),mh,
     .                wsp(ifree),lfree, iwsp, iexph, ns, iflag )
         iexph = ifree + iexph - 1
         nscale = nscale + ns
      else
*---     uniform rational Chebyshev approximation ...
         iexph = ifree
         do i = 1,mx
            wsp(iexph+i-1) = ZERO
         enddo
         wsp(iexph) = ONE
         call ZNCHBV(mx,sgn*t_step,wsp(ih),mh,wsp(iexph),wsp(ifree+mx))
      endif
 402  continue
* 
*---  error estimate ...
* 
      if ( k1.eq.0 ) then
         err_loc = tol
      else
         p1 = ABS( wsp(iexph+m) )   * beta
         p2 = ABS( wsp(iexph+m+1) ) * beta * avnorm
         if ( p1.gt.10.0d0*p2 ) then
            err_loc = p2
            xm = 1.0d0/DBLE( m )
         elseif ( p1.gt.p2 ) then
            err_loc = (p1*p2)/(p1-p2)
            xm = 1.0d0/DBLE( m )
         else
            err_loc = p1
            xm = 1.0d0/DBLE( m-1 )
         endif
      endif
*
*---  reject the step-size if the error is not acceptable ...
*   
      if ( (k1.ne.0) .and. (err_loc.gt.delta*t_step*tol) .and.
     .     (mxreject.eq.0 .or. ireject.lt.mxreject) ) then
         t_old = t_step
         t_step = gamma * t_step * (t_step*tol/err_loc)**xm
         p1 = 10.0d0**(NINT( LOG10( t_step )-SQR1 )-1)
         t_step = AINT( t_step/p1 + 0.55d0 ) * p1
         if ( itrace.ne.0 ) then
            print*,'t_step =',t_old
            print*,'err_loc =',err_loc
            print*,'err_required =',delta*t_old*tol
            print*,'stepsize rejected, stepping down to:',t_step
         endif
         ireject = ireject + 1
         nreject = nreject + 1
         if ( mxreject.ne.0 .and. ireject.gt.mxreject ) then
            print*,"Failure in ZGEXPV: ---"
            print*,"The requested tolerance is too high."
            Print*,"Rerun with a smaller value."
            iflag = 2
            return
         endif
         goto 401
      endif
*
*---  now update w = beta*V*exp(t_step*H)*e1 and the hump ...
*
      mx = mbrkdwn + MAX( 0,k1-1 )
      hij = DCMPLX( beta )
      
      CALL cublas_zdscal( ham_dim, ZERO, dev_w, 1)
      DO i = 1, mx
        CALL cublas_zaxpy(ham_dim,wsp(iexph+i-1),dev_v(i),1,dev_w,1)
      END DO
      CALL cublas_zdscal(n, hij, dev_w, 1)
      beta = cublas_dznrm2(n, dev_w, 1)
      
      hump = MAX( hump, beta )
*
*---  suggested value for the next stepsize ...
*
      t_new = gamma * t_step * (t_step*tol/err_loc)**xm
      p1 = 10.0d0**(NINT( LOG10( t_new )-SQR1 )-1)
      t_new = AINT( t_new/p1 + 0.55d0 ) * p1

      err_loc = MAX( err_loc,rndoff )
*
*---  update the time covered ...
*
      t_now = t_now + t_step
*
*---  display and keep some information ...
*
      if ( itrace.ne.0 ) then
         print*,'integration',nstep,'---------------------------------'
         print*,'scale-square =',ns
         print*,'step_size =',t_step
         print*,'err_loc   =',err_loc
         print*,'next_step =',t_new
      endif

      step_min = MIN( step_min, t_step )
      step_max = MAX( step_max, t_step )
      s_error = s_error + err_loc
      x_error = MAX( x_error, err_loc )
      if ( mxstep.eq.0 .or. nstep.lt.mxstep ) goto 100
      iflag = 1

 500  continue

      ! Here w = exp(a * ta(tt)) * v;
      cuda_status = cuda_memcpy_c2fort_real(w_seq((tt-1)*n+1), dev_w, 
     .  n * szpx, 2)
      
      END DO ta_output
      
      ! cuda, deallocate memory
      CALL cuda_free(dev_w)
      DO i = 1, m + 2
        CALL cuda_free(dev_v(i))
      END DO
      CALL hamvec_cuda3_term( cublas_handle )
      
      iwsp(1) = nmult
      iwsp(2) = nexph
      iwsp(3) = nscale
      iwsp(4) = nstep
      iwsp(5) = nreject
      iwsp(6) = ibrkflag
      iwsp(7) = mbrkdwn

      wsp(1)  = DCMPLX( step_min )
      wsp(2)  = DCMPLX( step_max )
      wsp(3)  = DCMPLX( 0.0d0 )
      wsp(4)  = DCMPLX( 0.0d0 )
      wsp(5)  = DCMPLX( x_error )
      wsp(6)  = DCMPLX( s_error )
      wsp(7)  = DCMPLX( tbrkdwn )
      wsp(8)  = DCMPLX( sgn*t_now )
      wsp(9)  = DCMPLX( hump/vnorm )
      wsp(10) = DCMPLX( beta/vnorm )
      
      END SUBROUTINE HAMVEC_ZGEXPV_W_CUDA_PROFILE
*----------------------------------------------------------------------|

