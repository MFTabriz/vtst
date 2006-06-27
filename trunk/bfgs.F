#include "symbol.inc"
!**********************************************************************
!
! Module which implements our full bfgs method
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE bfgs
    USE prec
    USE lattice

    IMPLICIT NONE
    private
    public :: bfgs_step, bfgs_init  !call bfgs_init from opt_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step ! ...
    LOGICAL :: fdstep, optflag_local

!**********************************************************************
!
! full bfgs method
!
!**********************************************************************

    subroutine bfgs_init()
    ! load the variables: dR, maxmove, memory
    ! initialize flag: finite_diff_step_flag
    end subroutine bfgs_init

    subroutine bfgs_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: posion_old,force_old
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,avg_force
      REAL(q) :: finite_step,step_size

      if (finite_difference_step_flag) then
        finite_difference_step_flag=.false.
        ! check for reset of direction
        a1 = abs(sum(force,force_old))
        a2 = sum(force_old,force_old)
        if (a1.lt.0.5_q*a2) then
           ! find new direct
        else
           ! the hessian to Idenity
        endif

        posion_old = posion
        force_old = force
        ! compute Ho*g
!**********************************************************************
!
! full bfgs method
!
!**********************************************************************

    subroutine bfgs_init()
    ! load the variables: dR, maxmove, memory
    ! initialize flag: finite_diff_step_flag
    end subroutine bfgs_init

    subroutine bfgs_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: posion_old,force_old
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,avg_force
      REAL(q) :: finite_step,step_size

      if (finite_difference_step_flag) then
        finite_difference_step_flag=.false.
        ! check for reset of direction
        a1 = abs(sum(force,force_old))
        a2 = sum(force_old,force_old)
        if (a1.lt.0.5_q*a2) then
           ! find new direct
        else
           ! the hessian to Idenity
        endif

        posion_old = posion
        force_old = force
        ! compute Ho*g

        direction = return_unit(direction)

        ! finite step down force
        posion = posion+direction*step_size
      else
        finite_difference_step_flag=.true.
        ! calculate curvature down direction
        fp1 = sum(force_old*direction)
        fp2 = sum(force*direction)
        avg_force = (fp1+fp2)/2.0_q
        curvature = (fp1-fp2)/finite_step
        if (curvature<0.0_q) then
          step_size=maxmove
        else
          step_size = avg_force/curvature
          if(abs(step_size)>maxmove) then
            step_size=sign(maxmove,step_size)
          else
            step_size = step_size-(finite_step/2.0)
          endif
        endif
        posion = posion+(direction*step_size)
        ! update the Hessian

      endif

    return


    end subroutine bfgs_step


