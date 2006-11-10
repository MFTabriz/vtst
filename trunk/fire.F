#include "symbol.inc"
!**********************************************************************
!
! Module for Fast Inertial Relaxation Engine (FIRE)
! from Erik Bitzek,PHYSICAL REVIEW LETTER  97, 170201 (2006)
!
! Version 1.01, Nov 2006
!
!**********************************************************************

  MODULE fire
    USE prec
    USE lattice

    IMPLICIT NONE
    SAVE
    private
    public :: fire_step, fire_init  !call fire_init from opt_init

    INTEGER :: nions,iu6,Nsteps,Nmin
    REAL(q),ALLOCATABLE :: step(:,:),velocity(:,:),R(:,:)
    REAL(q) :: maxmove,dt,dtmax
    REAL(q) :: finc,fdec,fa,alpha,alpha_start


!**********************************************************************
!
! Fast Inertial Relaxation Engine (FIRE) method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE fire_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      REAL(q) :: Power
      LOGICAL :: optflag

      R = posion
      
      optflag = .false. ! gives control back to the method

      ! Convert the position into Cartesian coordinates
      call dirkar(nions,R,latt_a)


      Power = SUM(force*velocity)
      IF (SUM(force*velocity) .GT. 0.0_q) THEN
        velocity = (1.0_q-alpha)*velocity+ &
          alpha*SUM(velocity*force)*force/SUM(force*force)
        IF (Nsteps .GT. Nmin) THEN
          ! Incease time step and decrease alpha
          dt = MIN(dt*finc,dtmax)
          alpha = alpha*fa
        ENDIF
        Nsteps = Nsteps+1
      ELSE
        ! reset alpha, velocity and decrease dt
        velocity = velocity*0.0_q
        alpha = alpha_start
        dt = dt*fdec
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

    END SUBROUTINE fire_step

!**********************************************************************
!
! Fire init
!
!**********************************************************************

    SUBROUTINE fire_init(T_INFO,IO)
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

      maxmove = 0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''MAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      dt = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'TIMESTEP','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''TIMESTEP'' from file INCAR.'
         STOP
      ENDIF
      
      dtmax = 1.0_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'TIMESTEPMAX','=','#',';','F', &
     &            IDUM,dtmax,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''TIMESTEPMAX'' from file INCAR.'
         STOP
      ENDIF
      
      finc = 1.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'TIMEINC','=','#',';','F', &
     &            IDUM,finc,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''TIMEINC'' from file INCAR.'
         STOP
      ENDIF
      
      fdec = 0.5_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'TIMEDEC','=','#',';','F', &
     &            IDUM,fdec,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''TIMEDEC'' from file INCAR.'
         STOP
      ENDIF
      
      fa = 0.99_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ALPHADEC','=','#',';','F', &
     &            IDUM,fa,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''ALPHADEC'' from file INCAR.'
         STOP
      ENDIF
      
      alpha_start = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ALPHASTART','=','#',';','F', &
     &            IDUM,alpha_start,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''ALPHASTART'' from file INCAR.'
         STOP
      ENDIF

      Nmin = 5
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'INMIN','=','#',';','I', &
     &            Nmin,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''INMIN'' from file INCAR.'
         STOP
      ENDIF
      
      
  
      ! initlize vectors
      ALLOCATE(step(3,nions),velocity(3,nions),R(3,nions))

      velocity = 0._q
      Nsteps = 0

    END SUBROUTINE fire_init

  END MODULE fire
