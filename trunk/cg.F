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
    REAL(q),ALLOCATABLE :: Fold(:,:),R(:,:)
    REAL(q),ALLOCATABLE :: direction(:,:),dir_unit(:,:)

    REAL(q) :: finite_step,maxmove
    LOGICAL :: fdstep

!**********************************************************************
!
! Conjugate gradient method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE cg_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      LOGICAL optflag
      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)

      REAL(q) :: a1,a2,gam,step_size
      REAL(q) :: fp1,fp2,Favg,curvature

      R = posion

      IF (fdstep) THEN
        fdstep = .false.
        optflag = .true.
        a1 = ABS(SUM(force*Fold))
        a2 = SUM(Fold*Fold)
        IF ((a1 .LT. 0.5_q*a2) .AND. (a2 .NE. 0.0_q)) THEN
          gam = SUM(force*(force-Fold))/a2
        ELSE
          gam = 0.0_q
        END IF
        direction = force+direction*gam
        dir_unit = direction/SQRT(SUM(direction*direction))
        Fold = force
        ! Move from the original configuration
        R = R+dir_unit*finite_step
      ELSE
        fdstep = .true.
        optflag = .false.
        fp1 = SUM(Fold*dir_unit)
        fp2 = SUM(force*dir_unit)
        curvature = (fp1-fp2)/finite_step
        IF (curvature .LT.  0.0_q) THEN
         step_size = maxmove
        ELSE
          Favg = 0.5_q*(fp1+fp2)
          step_size = Favg/curvature
          IF (ABS(step_size) .GT. maxmove) THEN
            step_size = SIGN(maxmove,step_size)
          ELSE
            step_size = step_size-0.5_q*finite_step   ! (*)
          END IF
        END IF
! Move now from the configuration after the fd_step, so (*) has a "-" sign
        R = R+dir_unit*step_size
      END IF

      ! convert position back to direct cord
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

    END SUBROUTINE cg_step
!**********************************************************************
!
! Conjugate gradient initilizer
!
!**********************************************************************

    SUBROUTINE cg_init(T_INFO,IO)
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
      iu6=IO%IU6

      ! read in varables used for cg

      finite_step=0.001_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'RSTEP','=','#',';','F', &
     &            IDUM,finite_step,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''RSTEP'' from file INCAR.'
         STOP
      ENDIF

      maxmove=0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''MAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      ! initialize the variables
 
      ALLOCATE(Fold(3,nions),R(3,nions))
      ALLOCATE(direction(3,nions),dir_unit(3,nions))

      fdstep=.true.  
    END SUBROUTINE cg_init

  END MODULE cg
