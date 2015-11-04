!=======================================================================
! tictac( clicker, sys_tic, cpu_tic, counter, timer )
!
! record the elapsed time
!
! Input:
!   clicker   (integer)   status of the switcher
!   sys_tic   (integer)   tic time of the system counter
!   cpu_tic   (double )   tic time of the cpu time
!
! Output:
!   counter   (integer)   number of the tictacs
!   timer(1)  (double )   accumulative cpu time
!   timer(2)  (double )   accumulative wall time
!=======================================================================
! Example:
!  INTEGER, PARAMETER  ::  ttn = 10  ! number of profiles
!  INTEGER             ::  tti       ! index  of profiles
!  INTEGER             ::  tt_click(ttn)   = 0, tt_sys(ttn)
!  REAL*8              ::  tt_cpu(ttn)
!  INTEGER             ::  tt_counter(ttn) = 0
!  REAL*8              ::  tt_timer(2*ttn) = 0.0D0
!  ! ******
!  ! Start of tictac
!  ! ******
!  tti = 1
!  CALL tictac(tt_click(tti), tt_sys(tti), tt_cpu(tti), &
!  & tt_counter(tti), tt_timer(2*tti-1) )
!  ! ******
!  ! Code to profile
!  ! ******
!  ! ******
!  !  End  of tictac, duplicated code as the start
!  ! ******
!  tti = 1
!  CALL tictac(tt_click(tti), tt_sys(tti), tt_cpu(tti), &
!  & tt_counter(tti), tt_timer(2*tti-1) )
!  ! Print the profile
!  WRITE(*,*) '=================================================='
!  WRITE(*,*) 'tictac profile for ******:'
!  ! f90 write
!  WRITE(*,*) '      index     counter             wal_time&
!  &                  cpu_time'
!  ! f77 write
!      WRITE(*,*) '          i     counter             wal_time
!     .        cpu_time'
!  DO i = 1, ttn
!    WRITE(*,*) i, tt_counter(i), tt_timer(2*i-1), tt_timer(2*i)
!  END DO
!  WRITE(*,*) '=================================================='
!=======================================================================
SUBROUTINE tictac( clicker, sys_tic, cpu_tic, counter, timer )

  INTEGER   clicker, sys_tic, counter
  REAL*8    cpu_tic, timer(2)
  
  INTEGER   sys_tac, sys_rate
  REAL*8    cpu_tac
  
  SELECT CASE (clicker)
    CASE (0)
      ! switch on
      clicker = 1
      CALL SYSTEM_CLOCK( COUNT = sys_tic )
      CALL CPU_TIME( cpu_tic )
    CASE (1)
      ! switch off
      CALL SYSTEM_CLOCK( COUNT = sys_tac, COUNT_RATE = sys_rate )
      CALL CPU_TIME( cpu_tac )
      
      counter  = counter + 1
      timer(1) = timer(1) + ( DBLE(sys_tac - sys_tic) / DBLE(sys_rate) )
      timer(2) = timer(2) + ( cpu_tac - cpu_tic )
      
      clicker = 0
  END SELECT

END SUBROUTINE tictac

