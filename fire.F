#include "symbol.inc"
!**********************************************************************
!
! Module for Fast Inertial Relaxation Engine (FIRE)
! from Erik Bitzek, Phys. Rev. Lett. 97, 170201 (2006)
!
!**********************************************************************

  MODULE fire
    USE prec
    USE lattice

    IMPLICIT NONE
    SAVE
    private
    public :: fire_step, fire_init  !call fire_init from opt_init

    INTEGER :: nions, iu0, iu6, Nsteps, Nmin, ISIF_local
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: step, velocity, R
    REAL(q) :: cvelocity(3,3), cstep(3,3)
    REAL(q) :: maxmove, dt, dtmax
    REAL(q) :: finc, fdec, fadec, alpha, alpha_start
    LOGICAL :: qmflag, cell_flag

!**********************************************************************
! Fast Inertial Relaxation Engine (FIRE) method
!**********************************************************************

  CONTAINS
    SUBROUTINE fire_step(optflag, posion, toten, force, hstress, latt_a, latt_b)

      REAL(q) :: posion(3,nions), toten, force(3,nions)
      REAL(q) :: hstress(3,3), latt_a(3,3), latt_b(3,3)
      REAL(q) :: Power
      LOGICAL :: optflag

      R = posion
      
      optflag = .false. ! give control back to the method

      ! Convert the position into Cartesian coordinates
      call dirkar(nions,R,latt_a)

      If (SUM(velocity) .NE. 0.0_q)THEN
        Power = SUM(force*velocity)
        IF (cell_flag) Power = SUM(force*velocity)+SUM(hstress*cvelocity)
        IF (Power .GT. 0.0_q) THEN

          IF(qmflag) THEN
            ! qm step -> keep velocity in direction of force
            velocity = SUM(velocity*force)*force/SUM(force*force)
            IF (cell_flag) cvelocity = SUM(cvelocity*hstress)*hstress/SUM(hstress*hstress)
          ELSE
            ! fire -> keeps most of the intial velocity
            velocity = (1.0_q-alpha)*velocity+ &
              alpha*force/SQRT(SUM(force*force))*SQRT(SUM(velocity*velocity))
            IF (cell_flag) THEN
              cvelocity = (1.0_q-alpha)*cvelocity+ &
                alpha*hstress/SQRT(SUM(hstress*hstress))*SQRT(SUM(cvelocity*cvelocity))
            ENDIF
          ENDIF

          IF (Nsteps .GE. Nmin) THEN
            ! Increase time step and decrease alpha
            dt = MIN(dt*finc,dtmax)
            alpha = alpha*fadec
          ENDIF
          Nsteps = Nsteps+1
        ELSE
          ! reset alpha, velocity and decrease dt
          velocity = velocity*0.0_q
          IF (cell_flag) cvelocity = cvelocity*0.0_q
          alpha = alpha_start
          dt = dt*fdec
          Nsteps = 0
        ENDIF
      ENDIF

      ! Euler step
      velocity = velocity+(dt*force)
      step = dt*velocity
      IF (SQRT(SUM(step*step)) .GE. maxmove) THEN
        step = maxmove*step/SQRT(SUM(step*step))
      ENDIF

      R = R+step
      ! convert position back to direct coordinates
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

      IF (cell_flag .OR. (ISIF_local==3)) THEN
        cvelocity = cvelocity+(dt*hstress)
        cstep = dt*cvelocity
        IF (SQRT(SUM(cstep*cstep)) .GE. maxmove) THEN
          cstep = maxmove*cstep/SQRT(SUM(cstep*cstep))
        ENDIF
        latt_a = latt_a + cstep
      ENDIF

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, TIMESTEP ',dt
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, alpha    ',alpha
      END IF



    END SUBROUTINE fire_step

!**********************************************************************
! Fire init
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
      
      dtmax = 1.0_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FTIMEMAX','=','#',';','F', &
     &            IDUM,dtmax,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FTIMEMAX'' from file INCAR.'
         STOP
      ENDIF
      
      finc = 1.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FTIMEINC','=','#',';','F', &
     &            IDUM,finc,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FTIMEINC'' from file INCAR.'
         STOP
      ENDIF
      
      fdec = 0.5_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FTIMEDEC','=','#',';','F', &
     &            IDUM,fdec,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FTIMEDEC'' from file INCAR.'
         STOP
      ENDIF
      
      fadec = 0.99_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FALPHADEC','=','#',';','F', &
     &            IDUM,fadec,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FALPHADEC'' from file INCAR.'
         STOP
      ENDIF
      
      alpha_start = 0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FALPHA','=','#',';','F', &
     &            IDUM,alpha_start,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FALPHA'' from file INCAR.'
         STOP
      ENDIF

      Nmin = 5
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FNMIN','=','#',';','I', &
     &            Nmin,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FNMIN'' from file INCAR.'
         STOP
      ENDIF

      qmflag = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'QMFLAG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,qmflag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''QMFLAG'' from file INCAR.'
         STOP
      ENDIF

      ! initlize vectors
      ALLOCATE(step(3,nions),velocity(3,nions),R(3,nions))

      velocity = 0._q
      cvelocity = 0._q
      Nsteps = 0
      alpha = alpha_start

      ! determines if the cell should change
      cell_flag=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LNEBCELL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,cell_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LNEBCELL'' from file INCAR.'
        STOP
      ENDIF

      ISIF_local = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ISIF','=','#',';','I', &
     &            ISIF_local,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ISIF'' from file INCAR.'
        STOP
      ENDIF

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, MAXMOVE  ',maxmove
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, TIMESTEP ',dt
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, FTIMEMAX ',dtmax
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, FTIMEINC ',finc
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, FTIMEDEC ',fdec
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, FALPHADEC',fadec
        WRITE(iu6,'(A5,A,F12.5)') 'OPT:',' FIRE, FALPHA   ',alpha_start
        WRITE(iu6,'(A5,A,I6)')    'OPT:',' FIRE, FNMIN    ',Nmin
      END IF

    END SUBROUTINE fire_init

  END MODULE fire
