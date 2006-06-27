#include "symbol.inc"
!**********************************************************************
!
! Module for Steepest dencent
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE sd
    USE prec
    USE lattice

    IMPLICIT NONE
    private 
    public :: sd_step, sd_init  !call sd_init from opt_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    LOGICAL :: fdstep, optflag_local


!**********************************************************************
!
! Steepest-Decent method
!
!**********************************************************************

    subroutine sd_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: step

      step = force * alpha
      if (sqrt(sum(step*step)).ge.maxmove then
        step = maxmove*return_unit(step)
      endif

      posion = posion+step

    return ! do i need this?
    end subroutine sd_step


!**********************************************************************
!
! Steepest-Decent init
!
!**********************************************************************

    subroutine sd_init()

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM
      REAL(q),SAVE alpha,maxmove
  
      ! loading variables for sd step

      alpha=0.1_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'FALPHA','=','#',';','F', &
     &            IDUM,alpha,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''alpha'' from file INCAR.'
         STOP
      ENDIF

      maxmove=0.2_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF


    end subroutine sd_init()

