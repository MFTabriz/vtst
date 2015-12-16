#include "symbol.inc"
!**********************************************************************
!
! Module for Quick-Min 
!
!**********************************************************************

  MODULE qm
    USE prec
    USE lattice

    IMPLICIT NONE
    SAVE
    private
    public :: qm_step, qm_init  !call qm_init from opt_init

    INTEGER :: nions, iu0, iu6, ISIF_local
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: step, velocity, R
    REAL(q) :: cvelocity(3,3), cstep(3,3)
    REAL(q) :: maxmove, dt
    LOGICAL :: cell_flag

  CONTAINS

!**********************************************************************
! Quick-Min method
!**********************************************************************

    SUBROUTINE qm_step(optflag,posion,toten,force,hstress,latt_a,latt_b)

      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: hstress(3,3),latt_a(3,3),latt_b(3,3)
      LOGICAL :: optflag

      R = posion
      
      optflag = .false. ! gives control back to the method

      ! Convert the position into Cartesian coordinates
      call dirkar(nions,R,latt_a)

      ! remove antiparallel components of the velocity along the force
!      WHERE (velocity*force<0) velocity=0._q
      ! decided to take this out because it is based upon the coordinate system

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
      ! convert position back to direct coordinates
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

      IF (cell_flag .OR. (ISIF_local==3)) THEN
        IF ((SUM(hstress*cvelocity)) .GE. 0.0_q) THEN
          cvelocity = SUM(cvelocity*hstress)*hstress/SUM(hstress*hstress)
        ELSE
          cvelocity = cvelocity*0.0_q
        ENDIF

        ! Euler step
        cvelocity = cvelocity+(dt*hstress)
        cstep = dt*cvelocity
        IF (SQRT(SUM(cstep*cstep)) .GE. maxmove) THEN
          cstep = maxmove*cstep/SQRT(SUM(cstep*cstep))
        ENDIF

        latt_a = latt_a + cstep
      ENDIF

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
      cvelocity=0._q

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' QM, MAXMOVE ',maxmove
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' QM, TIMESTEP',dt
      END IF

      ! determines of cell should change in NEB
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

    END SUBROUTINE qm_init

  END MODULE qm
