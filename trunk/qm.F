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
    REAL(q),ALLOCATABLE :: step(:,:),velocity(:,:)
    REAL(q) :: maxmove,dt


!**********************************************************************
!
! Qucik-min method
!
!**********************************************************************

    SUBROUTINE qm_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL :: optflag
      
      optflag = .false. ! gives control back to the method

      IF ((SUM(force*velocity)) .GE. 0.0_q) THEN
        velocity = SUM(velocity*force)*force/SQRT(SUM(force*force))
      ELSE
        velocity = velocity*0.0_q
      ENDIF

      ! Euler step
      velocity = velocity+(dt*force)
      step = dt*velocity
      IF (SQRT(SUM(step*step)) .GE. maxmove) THEN
        step = maxmove*step/SQRT(SUM(step*step))
      ENDIF

      posion = posion+step

    RETURN 
    END SUBROUTINE qm_step

!**********************************************************************
!
! Qucik-min init
!
!**********************************************************************

    SUBROUTINE qm_init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions=T_INFO%nions
      iu6=IO%IU6

      maxmove = 0.2_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,RMAXMOVE,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      dt = 0.01_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'TIMESTEP','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''TIMESTEP'' from file INCAR.'
         STOP
      ENDIF
  
      ! initlize vectors
      ALLOCATE(step(3,nions),velocity(3,nions))

    END SUBROUTINE qm_init

  END MODULE qm
