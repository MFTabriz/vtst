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
    REAL(q),ALLOCATABLE :: direction(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)
    REAL(q),ALLOCATABLE :: invH(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step,maxmove
    LOGICAL :: fdstep,dfp

!**********************************************************************
!
! full bfgs method
!
!**********************************************************************

    SUBROUTINE bfgs_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      REAL(q) :: a1,a2,curvature
      REAL(q) :: fp1,fp2,Favg,step_size

      nions=tinfo_nions
      iu6=tinfo_iu6
      dir2car=latt_a
      car2dir=latt_b

      IF (fdstep) THEN
        fdstep = .false.
        optflag = .true. ! keeps contol so it can do finite step
        ! check for reset of direction
        a1 = ABS(SUM(F*Fold))
        a2 = SUM(Fold*Fold)
        IF (a1 .LT. 0.5_q*a2) THEN
           ! find new direct
           !updateing the hessian
           CALL invH_update() 
        ELSE
           ! reset down the force
           ! reset the hessian
           !invH = idenity
        ENDIF

        Rold = posion
        Fold = force
        ! compute Ho*g


        ! direction = Ho*direction !Ho is idenity matrix
        direction = return_unit(direction)
        ! finite step down force
        posion = posion + direction*finite_step
      ELSE
        fdstep = .true.
        optflag = .false. ! gives control back to method
        ! calculate curvature down direction
        fp1 = SUM(Fold*direction)
        fp2 = SUM(force*direction)
        curvature = (fp1-fp2)/finite_step
        IF (curvature .LT. 0.0_q) THEN
          step_size=maxmove
        ELSE
          Favg = 0.5_q*(fp1+fp2)
          step_size = Favg/curvature
          IF(ABS(step_size) .GT. maxmove) THEN
            step_size=sign(maxmove,step_size)
          ELSE
            step_size = step_size-0.5_q*finite_step  !(*)
          ENDIF
        ENDIF
! Move now from the configuration after the fd_step, so (*) has a "-" sign
        posion = posion+(direction*step_size)
      ENDIF
    RETURN
    END SUBROUTINE bfgs_step

!**********************************************************************
! Update the Inv HESSIAN
! William Press's Numerical Recipes in C p 427
!**********************************************************************
 ! need to get function and check on reshaping
    SUBROUTINE invH_update()
      REAL(q) :: a,b
          a = SUM(change_in_R*change_in_F)
          W = outter_product(change_in_R,change_in_R) !matrix
          A = matrixmultiply(invH,change_in_F) !vector
          AA = outerproduct(A,A)  !matrix
          b = vdot(change_in_F,A)  !scaler
          ! DFP method
          invH = invH+W/a-AA/b
          IF(.NOT. DFP) THEN
            u = change_in_R/a-A/b  !vector
            uu = outerproduct(u,u) !matrix
            invH = invH+W/a-AA/b+b*uu
          ENDIF
    END SUBROUTINE invH_update

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
      ALLOCATE(Fold(3,nions),Rold(3,nions),direction(3,nions))
      ALLOCATE(invH(3*nions,3*nions))

     ! set invH as indenty
      fdstep = .true.
    END SUBROUTINE bfgs_init
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
      SUBROUTINE set_pbc(v1)
      REAL(q) :: v1(3,nions)
      CALL kardir(nions,v1,car2dir)
      v1=MOD(v1+100.5_q,1._q)-0.5_q
      CALL dirkar(nions,v1,dir2car)
      END SUBROUTINE set_pbc
!======================================================================
! Returns a unit vector along v1
!======================================================================
      FUNCTION return_unit(V1)
      real(q) :: v1(3,nions)
      real(q),dimension(3,nions) :: return_unit
      return_unit=v1*(1._q/SQRT(SUM(v1*v1)))
      END FUNCTION return_unit
!======================================================================
! Sets V1 to be a unit vector
!======================================================================
      SUBROUTINE set_unit(V1)
      REAL(q) :: v1(3,nions)
      v1=return_unit(v1)
      END SUBROUTINE set_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      FUNCTION vproj(v1,v2)
      REAL(q) :: v1(3,nions),v2(3,nions),proj(3,nions)
      proj=v2*SUM(v1*v2)/SQRT(SUM(v2*v2))
      END FUNCTION vproj

  END MODULE bfgs
