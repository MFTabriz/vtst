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

    REAL(q),ALLOCATABLE :: step(:,:),R(:,:)
    REAL(q) :: alpha,maxmove


!**********************************************************************
!
! Steepest-Decent method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE sd_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag
      
      optflag = .false. ! gives control back to the method

      R = posion

      ! Convert the position into Cartesian coordinates
      call dirkar(nions,R,latt_a)

      step = force * alpha
      IF (SQRT(SUM(step*step)) .GT. maxmove) THEN
        step = maxmove*step/SQRT(SUM(step*step))
      ENDIF

      R = R+step

      ! convert position back to direct cord
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

    END SUBROUTINE sd_step


!**********************************************************************
!
! Steepest-Decent init
!
!**********************************************************************

    SUBROUTINE sd_init(T_INFO,IO)
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
      iu6=IO%IU6 

      ! loading variables for sd step

      alpha = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FALPHA','=','#',';','F', &
     &            IDUM,alpha,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''alpha'' from file INCAR.'
         STOP
      ENDIF

      maxmove = 0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      ! initilize the vectors
      ALLOCATE(step(3,nions),R(3,nions))

    END SUBROUTINE sd_init

  END MODULE sd
