#include "symbol.inc"
!**********************************************************************
!
! Module which implements our optimizers
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE cg
    USE prec
    USE lattice

    IMPLICIT NONE
    private 
    public :: cg_step, cg_init  !call cg_init from opt_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step ! ...
    LOGICAL :: fdstep, optflag_local

!**********************************************************************
!
! Conjugate gradient method
!
!**********************************************************************

    subroutine cg_init()
    ! load the variables: dR, maxmove
    ! initialize flag: finite_diff_step_flag
    end subroutine cg_init

    subroutine cg_step(posion,toten,force)
!      LOGICAL,INTENT(INOUT) :: fd_step
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q),PARAMETER :: gamma=0.5_q
      REAL(q) :: c1,c2,f1d,f2d,fr,cr,step,s

!      IF (fd_step) THEN
      IF (finite_difference_step_flag) THEN
        finite_difference_step_flag=.false.
        f1=Feff
        c1=ABS(SUM(f1*fco))
        c2=SUM(fco*fco)
        IF (c1 < gamma*c2) THEN
          s=SUM(f1*(f1-fco))/c2
        ELSE
          s=0.0_q
        END IF
        dc=f1+s*dco
        dcu=dc/SQRT(SUM(dc*dc))
        fco=f1
        dco=dc
! Move from the original configuration
        R=R0+dcu*dR
      ELSE
!        fd_step=.true.
        finite_difference_step_flag=.true.
        f2=Feff
        f1d=SUM(F1*dcu)
        f2d=SUM(f2*dcu)
        cr=(f1d-f2d)/dR
        IF (cr < 0.0_q) THEN
         step=maxmove
        ELSE
          fr=0.5_q*(f1d+f2d)
          step=fr/cr
          IF (ABS(step) > maxmove) THEN
            step=SIGN(maxmove,step)
          ELSE
            step=step-0.5_q*dR     ! (*)
          END IF
        END IF
! Move now from the configuration after the fd_step, so (*) has a "-" sign
        R=R+dcu*step
      END IF

    RETURN


    end subroutine cg_step

