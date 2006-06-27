#include "symbol.inc"
!**********************************************************************
!
!
! Module for limited-memory bfgs
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE lbfgs
    USE prec
    USE lattice

    IMPLICIT NONE
    private 
    public :: lbfgs_step, lbfgs_init  !call opt_init from chain_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step ! ...
    LOGICAL :: fdstep, optflag_local

!**********************************************************************
!
! Limited-memory bfgs method
!
!**********************************************************************

    subroutine lbfgs_init()

      REAL(q),ALLOCATABLE,DIMENSION(:,:,:),SAVE :: change_in_force,change_in_posion

      INTEGER IDUM,IERR,N
      INTEGER,SAVE memory
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      memory=4
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'NLBFGSMEM','=','#',';','F', &
     &            memory,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''NBFGSMEM'' from file INCAR.'
         STOP
      ENDIF

      ! initialize the variables: vectors and matricies
      allocate(change_in_force(memory,3,Nions),change_in_posion(memory,3,Nions))

      fdstep=.TRUE.
    end subroutine lbfgs_init

    subroutine lbfgs_step()
      REAL(q) :: change_in_force,change_in_posion !matricies

      REAL(q) :: rho,q

      REAL(q) :: alpha,beta ! these are vectors
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,avg_force
      REAL(q) :: finite_step,step_size
      INTEGER :: itr,incr,bound,im,jm

      if (finite_difference_step_flag) then
        finite_difference_step_flag=.false.
        ! check for reset of direction
        a1 = abs(sum(force*force_old))
        a2 = sum(force_old*force_old)
        if (a1.lt.0.5_q*a2) then
           ! find new direct
           change_in_force(itr,:,:) = F-Fold
           change_in_posion(itr,:,:) = R-Rold  ! need boundary conditions
           rho(itr) = 1.0_q/sum(change_in_force(itr,:,:)*change_in_posion(itr,:,:))
           itr = itr+1
        else
           ! reset down the force
           itr = 1
           ! ? reset the change_in_pos/force vectors?
        endif

        posion_old = posion
        force_old = force
        ! compute Ho*g
        if (itr.lt.memory) then
          incr = 0
          bound = itr
        else
          incr = itr-memory
          bound = memory
        endif
        q = force
        do im=0,bound-1  ! change the range to 1..N
          jm = (bound-im)+incr
          alpha(bound-im) = rho(jm)*sum(change_in_posion(jm)*q)
          q = q-alpha(bound-im)*change_in_force(jm)
        enddo
        direction = Ho*q
        do im=0,bound-1
          jm = im+incr
          beta = rho(jm)*sum(change_in_force(jm)*direction
          direction = direction + change_in_posion*(alpha(im)-beta)
        enddo
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
      endif

    return

    end subroutine lbfgs_step


