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
    REAL(q) :: maxmove,dt,velocity


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

      REAL(q) :: step
      
      optflag = .false. ! gives control back to the method

      IF ((SUM(force*velocity)) .GE. 0.0_q) THEN
        velocity = vproj(velocity,force)
      ELSE
        velocity = velocity*0.0_q
      ENDIF

      ! Euler step
      velocity = velocity*dt*force!/mass
      step = dt*velocity
      IF (sqrt(sum(step*step)) .GE. maxmove) THEN
        step = maxmove*return_unit(step)
      ENDIF

      posion = posion+step

    RETURN 
    END SUBROUTINE qm_step

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
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF
  
  ! ls velocity
      velocity = 

    END SUBROUTINE qm_init
!**********************************************************************
!
! Vector Functions
!
!**********************************************************************

!======================================================================
! Sets a vector to have the smallest length consistent the the periodic
! boundary conditions. This should really be changed to use the
! Wigner-Sitz cell
!======================================================================
      subroutine set_pbc(v1)
      real(q) :: v1(3,nions)
      call kardir(nions,v1,car2dir)
      v1=mod(v1+100.5_q,1._q)-0.5_q
      call dirkar(nions,v1,dir2car)
      end subroutine set_pbc
!======================================================================
! Returns a unit vector along v1
!======================================================================
      function return_unit(V1)
      real(q) :: v1(3,nions)
      real(q),dimension(3,nions) :: return_unit
      return_unit=v1*(1._q/sqrt(sum(v1*v1)))
      end function return_unit
!======================================================================
! Sets V1 to be a unit vector
!======================================================================
      subroutine set_unit(V1)
      real(q) :: v1(3,nions)
      v1=return_unit(v1)
      end subroutine set_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      function vproj(v1,v2)
      real(q) :: v1(3,nions),v2(3,nions),proj(3,nions)
      proj=v2*sum(v1*v2)/sum(v2*v2)
      end function vproj

  END MODULE qm
