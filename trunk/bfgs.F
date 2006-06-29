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
    REAL(q) :: finite_step,maxmove
    LOGICAL :: fdstep,dfp

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

      F = force
      R = posion

      IF (fdstep) THEN
        fdstep = .false.
        optflag = .true. ! keeps contol so it can do finite step
        ! check for reset of direction
        a1 = ABS(SUM(F*Fold))
        a2 = SUM(Fold*Fold)
        IF (a1 .LT. 0.5_q*a2) THEN
           ! find new direct
           !updateing the hessian
        ELSE
           ! reset down the force
           ! ? reset the hessian
        ENDIF

        Rold = F
        Fold = F
        ! compute Ho*g
        IF (itr .LT. memory) THEN
          bound = itr
        ELSE
          bound = memory
        ENDIF
        direction = force
        DO im=1,bound
          jm = (bound+1-im)
          alpha(im) = rho(jm)*SUM(change_in_R(jm)*direction)
          direction = direction-(alpha(im)*change_in_F(jm))
        ENDDO
        ! direction = Ho*direction !Ho is idenity matrix
        DO im=1,bound
          beta = rho(im)*SUM(change_in_F(im)*direction
          direction = direction+change_in_R*(alpha(im)-beta)
        ENDDO
        direction = return_unit(direction)

        ! finite step down force
        R = R + direction*finite_step
      ELSE
        fdstep = .true.
        optflag = .false. ! gives control back to method
        ! calculate curvature down direction
        fp1 = SUM(Fold*direction)
        fp2 = SUM(F*direction)
        curvature = (fp1-fp2)/finite_step
        IF (curvaturR .LT. 0.0_q) THEN
          step_size=maxmove
        ELSE
          Favg = 0.5_q*(fp1+fp2)
          step_size = Favg/curvature
          IF(ABS(step_size) .GT. maxmove) THEN
            step_size=sign(maxmove,step_size)
          ELSE
            step_size = step_size-0.5_q*finite_step
          ENDIF
        ENDIF
        R = R+(direction*step_size)
      ENDIF

      posion(:,:) = R


    RETURN
    END SUBROUTINE bfgs_step

!**********************************************************************
!
! full bfgs init
!
!**********************************************************************

    SUBROUTINE bfgs_init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions = T_INFO%nions
      iu6 = IO%IU6

      ! read in the varables for bfgs

      finite_step = 0.001_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RSTEP','=','#',';','F', &
     &            IDUM,RSTEP,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RSTEP'' from file INCAR.'
         STOP
      ENDIF

      maxmove = 0.2_q
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''RMAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      dfp = .true.
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'RMAXMOVE','=','#',';','F', &
     &            IDUM,RDUM,CDUM,BFGSDFP,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''BFGSDFP'' from file INCAR.'
         STOP
      ENDIF


     ! intitilize vectors and matrix
    
      fdstep = .true.
    END SUBROUTINE bfgs_init
  END MODULE bfgs
