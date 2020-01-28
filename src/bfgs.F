#include "symbol.inc"
!**********************************************************************
!
! Module which implements the BFGS method
!
! To do 
!   - matrix operations
!   - boundary conditions
!
!**********************************************************************

  MODULE bfgs
    USE prec
    USE lattice

    IMPLICIT NONE
    SAVE
    PRIVATE
    PUBLIC :: bfgs_step, bfgs_init  ! CALL bfgs_init from opt_init

    INTEGER :: nions,iu6,itr,im,ij,ik
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: R,F,Rold,Fold,step
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: change_in_R,change_in_G,invH
    REAL(q),ALLOCATABLE,DIMENSION(:) :: direction

    ! temporary variables
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: SY,HyHy,UU
    REAL(q),ALLOCATABLE,DIMENSION(:) :: s,y,u,Hy

    REAL(q),DIMENSION(3,3) :: dir2car,car2dir
    REAL(q) :: finite_step,maxmove,init_curvature
    LOGICAL :: fdstep,dfp,line,reset_flag

!**********************************************************************
!
! BFGS method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE bfgs_step(optflag,posion,toten,force,latt_a,latt_b)

      REAL(q),DIMENSION(3,nions) :: posion,force
      REAL(q),DIMENSION(3,3) :: latt_a,latt_b
      REAL(q) :: toten
      LOGICAL optflag

      REAL(q) :: a1,a2,curvature
      REAL(q) :: fp1,fp2,Favg,step_size
      INTEGER :: im

      R = posion
      F = force
      dir2car=latt_a
      car2dir=latt_b

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A)') 'OPT:','INSIDE BFGS step function'
      ENDIF


      ! Convert the position into Cartesian coordinates
      CALL dirkar(nions,R,latt_a)

      IF (fdstep) THEN ! Finite difference step

        IF (iu6>=0) THEN
          WRITE(iu6,'(A5,A)') 'OPT:','taking fdstep'
        ENDIF

        fdstep = .FALSE.
        optflag = .TRUE. ! keep control to do finite difference step

        ! check for reset of direction
        a1 = ABS(SUM(F*Fold))
        a2 = SUM(Fold*Fold)
!        IF (a1 .LT. 0.5_q*a2) THEN
        IF ((a1 .GT. 0.5_q*a2) .OR. (a2 .EQ. 0.0_q)) THEN
           reset_flag=.TRUE.
           IF (iu6>=0) THEN
             WRITE(iu6,'(A5,A)') 'OPT:','reseting direction'
           ENDIF
        ENDIF

        ! reset the Hessian
        IF(reset_flag) THEN
          itr=0
          reset_flag=.FALSE.
          ! set invH identity matrix
          invH = 0.0_q
          DO im=1,3*nions
            invH(im,im)=1.0_q
          ENDDO
          invH=invH*init_curvature

        ELSE
          change_in_R = R-Rold
          ! apply periodic boundary conditions to these vector differences
          CALL set_pbc(change_in_R)
          change_in_G = F-Fold
          

          ! update the (inverse) Hessian
          CALL invH_update(change_in_R,change_in_G,invH) 
        ENDIF

        Rold = R
        Fold = force

        direction = RESHAPE(F,SHAPE=(/ 3*nions /))
        ! compute -Ho*g
        direction = MATMUL(invH,direction)
        step = RESHAPE(direction,SHAPE=(/3,nions/))

        IF(line) THEN
          ! line minimizer
          step = return_unit(step)
          ! finite step down force
          R = R + step*finite_step
        ELSE
          ! step_size based on Hessian
          IF(SQRT(ABS(SUM(step*step))) .GT. maxmove) THEN
            step = return_unit(step)*maxmove
          ENDIF
          R = R + step
          fdstep = .TRUE.
          optflag = .false. ! give control back to method
        ENDIF

      ! calculate Newton step based on curvature
      ELSE
        IF (iu6>=0) THEN
          WRITE(iu6,'(A5,A)') 'OPT:','taking big step'
        ENDIF
        fdstep = .TRUE.
        optflag = .false. ! give control back to method
        ! calculate curvature down direction
        fp1 = SUM(Fold*step)
        fp2 = SUM(F*step)
        curvature = (fp1-fp2)/finite_step
        IF (curvature .LT. 0.0_q) THEN

          IF (iu6>=0) THEN
            WRITE(iu6,'(A5,A)') 'OPT:','neg curv maxmove'
          ENDIF
          step_size=maxmove
        ELSE
          Favg = 0.5_q*(fp1+fp2)
          step_size = Favg/curvature
          IF(ABS(step_size) .GT. maxmove) THEN
            step_size=SIGN(maxmove,step_size)
          ELSE
            step_size = step_size-0.5_q*finite_step  !(*)
          ENDIF
        ENDIF
! Move now from the configuration after the fd_step, so (*) has a "-" sign
        R = R+(step*step_size)
        IF (iu6>=0) THEN
          WRITE(iu6,'(A5,A)') 'OPT:','taking bigstep done'
        ENDIF
      ENDIF

      ! convert position back to direct coordinates
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R
    END SUBROUTINE bfgs_step

!**********************************************************************
! Update the inverse Hessian
! William Press Numerical Recipes in C p 427
! TODO: need to get function and check on reshaping
!**********************************************************************

    SUBROUTINE invH_update(change_in_R,change_in_G,invH)
      REAL(q),DIMENSION(:,:) :: change_in_R,change_in_G,invH
!      REAL(q) :: s(:),y(:),u(:),Hy(:)
!      REAL(q) :: SY(:,:),HyHy(:,:),UU(:,:)

      REAL(q) :: a,b
      INTEGER :: im,jm,km
          

        IF (iu6>=0) THEN
          WRITE(iu6,'(A5,A)') 'OPT: ','IN HESS update'
        ENDIF
        s = RESHAPE(change_in_R,SHAPE=(/ 3*nions /))
        y = RESHAPE(change_in_G,SHAPE=(/ 3*nions /))
        a = SUM(s*y)           !scaler
        CALL outerproduct(s,s,SY)  !matrix
!        Hy = MATMUL(invH,y)     !vector
!        HyHy = outerproduct(Hy,Hy) !matrix
!        b = SUM(y*Hy)           !scaler
        ! DFP method
!        invH = invH+SY/a-HyHy/b
!        IF(.NOT. DFP) THEN
!          u = s/a-Hy/b  !vector
!          UU = outerproduct(u,u) !matrix
!          invH = invH+SY/a-HyHy/b+b*UU
!        ENDIF
        IF (iu6>=0) THEN
          WRITE(iu6,'(A5,A)') 'OPT: ','IN HESS update done'
        ENDIF
    END SUBROUTINE invH_update

!**********************************************************************
! BFGS init
!**********************************************************************

    SUBROUTINE bfgs_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,iu0
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions = T_INFO%nions
      iu0 = IO%IU0
      iu6 = IO%IU6

      ! read the variables for BFGS

      IF (iu6>=0) THEN
        WRITE(iu6,*) 'OPT: ','!!!!!!!!WARNING!!!!!!!!!!!! ' 
        WRITE(iu6,*) 'OPT: ','!!!!BFGS not complete!!!!!! '
        WRITE(iu6,*) 'OPT: ','Must change setting for IOPT'
      ENDIF
      STOP

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A)') 'OPT:','INSIDE BFGS init function'
      ENDIF


      finite_step = 0.005_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FDSTEP','=','#',';','F', &
     &            IDUM,finite_step,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''FDSTEP'' from file INCAR.'
         STOP
      ENDIF

      maxmove = 0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''MAXMOVE'' from file INCAR.'
         STOP
      ENDIF

      dfp = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'BFGSDFP','=','#',';','L', &
     &            IDUM,RDUM,CDUM,dfp,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''BFGSDFP'' from file INCAR.'
         STOP
      ENDIF

      init_curvature = 0.0001_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'BFGSINVCURV','=','#',';','F', &
     &            IDUM,init_curvature,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''BFGSINVCURV'' from file INCAR.'
         STOP
      ENDIF

      line = .TRUE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LINEMIN','=','#',';','L', &
     &            IDUM,RDUM,CDUM,line,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''LINEMIN'' from file INCAR.'
         STOP
      ENDIF

      ! intitilize vectors and matrix
      ALLOCATE(R(3,nions),F(3,nions))
      ALLOCATE(Fold(3,nions),Rold(3,nions))
      ALLOCATE(change_in_R(3,nions),change_in_G(3,nions))
      ALLOCATE(step(3,nions),direction(3*nions))
      ALLOCATE(invH(3*nions,3*nions))

      ! intitilize temporary variables
      ALLOCATE(SY(3*nions,3*nions),HyHy(3*nions,3*nions))
      ALLOCATE(UU(3*nions,3*nions),s(3*nions),y(3*nions))
      ALLOCATE(u(3*nions),Hy(3*nions))

      R = 0.0_q
      F = 0.0_q
      Rold = 0.0_q
      Fold = 0.0_q
      change_in_R = 0.0_q
      change_in_G = 0.0_q
      step = 0.0_q
      direction = 0.0_q

      ! set invH identity matrix
      invH = 0.0_q

      DO im=1,3*nions
        invH(im,im)=1.0_q
      ENDDO
      invH=invH*init_curvature

      fdstep = .TRUE.
      reset_flag = .TRUE.

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:','BFGS, FDSTEP',finite_step
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:','BFGS, MAXMOVE',maxmove
!        WRITE(iu6,'(A5,A,L14)') 'OPT:','BFGS, BFGSDFP',dfp
      END IF

      IF(iu6>=0)  WRITE(iu6,*) 'OPT:','done bfgs_init'

    END SUBROUTINE bfgs_init

!**********************************************************************
! Vector Functions
!**********************************************************************


!======================================================================
! Preforms outer product of two vectors
!======================================================================
      SUBROUTINE outerproduct(v1,v2,M1)
      INTEGER im,ik
      REAL(q) :: v1(3*nions),v2(3*nions)
      REAL(q),dimension(3*nions,3*nions) :: M1
      DO im = 1,3*nions
        DO ik = 1,3*nions
          M1(im,ik)=v1(im)*v2(ik)
        ENDDO
      ENDDO
      IF(iu6>=0)  WRITE(iu6,*) 'OPT:','in outterproduct'
      END SUBROUTINE outerproduct
!======================================================================
! Sets a vector to have the smallest length consistent with the periodic
! boundary conditions.
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
      REAL(q) :: v1(3,nions),v2(3,nions),vproj(3,nions)
      vproj=v2*SUM(v1*v2)/SUM(v2*v2)
      END FUNCTION vproj

  END MODULE bfgs
