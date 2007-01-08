#include "symbol.inc"
!**********************************************************************
!
! Module for Quick-min 
!
! Version 2.00, Jan 2007
!
!**********************************************************************

  MODULE qm
    USE prec
    USE lattice

    IMPLICIT NONE
    SAVE
    private
    public :: qm_step, qm_init  !call qm_init from opt_init

    INTEGER :: nions,iu0,iu6
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: step,velocity,R
    REAL(q) :: maxmove,dt

  CONTAINS

!**********************************************************************
! Quick-Min method
!**********************************************************************

    SUBROUTINE qm_step(optflag,posion,toten,force,latt_a,latt_b)

      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL :: optflag

      R = posion
      
      optflag = .false. ! gives control back to the method

      ! Convert the position into Cartesian coordinates
      call dirkar(nions,R,latt_a)

      ! remove antiparallel components of the velocity along the force
      WHERE (velocity*force<0) velocity=0._q

      IF ((SUM(force*velocity)) .GE. 0.0_q) THEN
        velocity = SUM(velocity*force)*force/SUM(force*force)
      ELSE
        velocity = velocity*0.0_q
      ENDIF

      ! Euler step
      velocity = velocity+(dt*force)
      step = dt*velocity
      IF (SQRT(SUM(step*step)) .GE. maxmove) THEN
        step = maxmove*step/SQRT(SUM(step*step))
      ENDIF

      R = R+step
      ! convert position back to direct cord
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

    END SUBROUTINE qm_step

!**********************************************************************
! Quick-Min init
!**********************************************************************

    SUBROUTINE qm_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions=T_INFO%nions
      iu0=IO%IU0
      iu6=IO%IU6

      maxmove = 0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''MAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      dt = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'TIMESTEP','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''TIMESTEP'' from file INCAR.'
         STOP
      ENDIF
  
      ! initialize vectors
      ALLOCATE(step(3,nions),velocity(3,nions),R(3,nions))
      velocity=0._q

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:','QM, MAXMOVE',maxmove
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:','QM, TIMESTEP',dt
      END IF

    END SUBROUTINE qm_init

  END MODULE qm
