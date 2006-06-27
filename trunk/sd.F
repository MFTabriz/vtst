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
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    LOGICAL :: fdstep, optflag_local


!**********************************************************************
!
! Steepest-Decent method
!
!**********************************************************************

    subroutine sd_init()
    ! load the variables: alpha, maxmove
    end subroutine sd_init()

    subroutine sd_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: step

      step = force * alpha
      if (sqrt(sum(step*step)).ge.maxmove then
        step = maxmove*return_unit(step)
      endif

      posion = posion+step

    return ! do i need this?
    end subroutine sd_step


