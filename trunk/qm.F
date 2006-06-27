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

    subroutine qm_init()
    ! load the variables: dt, maxmove
    ! set velocity
    end subroutine qm_init

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
