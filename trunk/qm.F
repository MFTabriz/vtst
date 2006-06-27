#include "symbol.inc"
!**********************************************************************
!
! Module for Quick-min 
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE qm
    USE prec
    USE lattice

    IMPLICIT NONE
    private
    public :: qm_step, qm_init  !call qm_init from opt_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    LOGICAL :: fdstep, optflag_local


!**********************************************************************
!
! Qucik-min method
!
!**********************************************************************

    subroutine qm_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: velocity,step

      if ((sum(force*velocity)).ge.0.0_q) then
        velocity = vproj(velocity,force)
      else
        velocity = velocity*0.0_q
      endif

      ! Euler step
      velocity = velocity*dt*force/mass
      step = dt*velocity
      if (sqrt(sum(step*step)).ge.maxmove) then
        step = maxmove*return_unit(step)
      endif

      posion = posion+step

    return ! do i need this?
    end subroutine qm_steP

!**********************************************************************
!
! Qucik-min init
!
!**********************************************************************

    subroutine qm_init()

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM
      REAL(q),SAVE maxmove,dt

      maxmove=0.2_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,RMAXMOVE,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      dt=0.01_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'TIMESTEP','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF
  
  ! set velocity

    end subroutine qm_init
