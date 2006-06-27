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

    subroutine bfgs_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: posion_old,force_old
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,avg_force
      REAL(q) :: finite_step,step_size

      if (fdstep) then
        fdstep=.false.
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
        fdstep=.true.
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

!**********************************************************************
!
! full bfgs init
!
!**********************************************************************

    subroutine bfgs_init()

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM
      REAL(q),SAVE finite_step,maxmove

      ! read in the varables for bfgs

      finite_step=0.001_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RSTEP','=','#',';','F', &
     &            IDUM,RSTEP,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RSTEP'' from file INCAR.'
         STOP
      ENDIF

      maxmove=0.2_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      fdstep = .TRUE.
    end subroutine bfgs_init
