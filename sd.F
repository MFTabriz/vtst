#include "symbol.inc"
!**********************************************************************
!
! Module for Steepest Descents
!
!**********************************************************************

  MODULE sd
    USE prec
    USE lattice

    IMPLICIT NONE
    PRIVATE 
    PUBLIC :: sd_step, sd_init  ! call sd_init from opt_init

    INTEGER :: nions,iu6

    REAL(q),ALLOCATABLE :: step(:,:),R(:,:)
    REAL(q) :: alpha,maxmove

!**********************************************************************
! Steepest-Descent method
!**********************************************************************

  CONTAINS
    SUBROUTINE sd_step(optflag,posion,toten,force,latt_a,latt_b)

      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag
      
      optflag = .false. ! give control back to the method

      R = posion

      ! Convert the position into Cartesian coordinates
      CALL dirkar(nions,R,latt_a)

      step = force*alpha
      IF (SQRT(SUM(step*step)) .GT. maxmove) THEN
        step = maxmove*step/SQRT(SUM(step*step))
      ENDIF
      R = R+step

      IF (iu6>=0) THEN
        WRITE(iu6,'(A4,A,F14.6)') 'OPT:','SD, step size',SQRT(SUM(step*step))
      END IF

      ! convert position back to direct coordinates
      CALL kardir(nions,R,latt_b)

      ! update posion
      posion = R

    END SUBROUTINE sd_step

!**********************************************************************
! Steepest-Descent init
!**********************************************************************

    SUBROUTINE sd_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,iu0
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions=T_INFO%nions
      iu0=IO%iu0
      iu6=IO%iu6

      ! loading variables for sd step

      alpha = 0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'SDALPHA','=','#',';','F', &
     &            IDUM,alpha,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''alpha'' from file INCAR.'
         STOP
      ENDIF

      maxmove = 0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''MAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:','SD, SDALPHA',alpha
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:','SD, MAXMOVE',maxmove
      END IF

      ! initialize the vectors
      ALLOCATE(step(3,nions),R(3,nions))

    END SUBROUTINE sd_init

  END MODULE sd
