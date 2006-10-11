#include "symbol.inc"
!**********************************************************************
!
! Module for limited-memory bfgs
!
! Version 1.01, June 2006
!
! TO DO
!   - PBC
!   - figure out loop exactly
!   - write shifting matrix
!
!**********************************************************************

  MODULE lbfgs
    USE prec
    USE lattice

    IMPLICIT NONE
    PRIVATE 
    PUBLIC :: lbfgs_step, lbfgs_init  !CALL opt_init from chain_init

    INTEGER :: nions,iu6,memory
    REAL(q),ALLOCATABLE ::R(:,:),Rold(:,:),Fold(:,:)
    REAL(q),ALLOCATABLE :: rho(:),alpha(:)
    REAL(q),ALLOCATABLE :: direction(:,:)
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:) :: change_in_F,change_in_R

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step,maxmove
    LOGICAL :: fdstep

!**********************************************************************
!
! Limited-memory bfgs method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE lbfgs_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      REAL(q) :: beta
      REAL(q) :: a1,a2,curvature
      REAL(q) :: fp1,fp2,Favg,step_size
      INTEGER :: itr,bound,im,jm

      R = posion

      ! Convert the position into Cartesian coordinates
      CALL dirkar(nions,R,latt_a)

      IF (fdstep) THEN
        fdstep = .FALSE.
        optflag = .TRUE. ! keeps contol so it can do finite step
        ! check for reset of direction
        a1 = ABS(SUM(force*Fold))
        a2 = SUM(Fold*Fold)
        IF (a1 .LT. 0.5_q*a2 .AND. a2 .ne. 0.0_q) THEN
           ! find new direct
           itr = itr+1
           IF (itr .LT. memory) THEN
             change_in_F(itr,:,:) = force-Fold
             change_in_R(itr,:,:) = R-Rold 
             ! apply periodic boundary conditions to these vector differences
             CALL set_pbc(change_in_R(itr,:,:))
             rho(itr) = 1.0_q/SUM(change_in_F(itr,:,:)*change_in_R(itr,:,:))
           ELSE
             DO im = 1,memory-1
               change_in_F(im,:,:) = change_in_F(im+1,:,:)
             ENDDO
             DO im = 1,memory-1
               change_in_R(im,:,:) = change_in_R(im+1,:,:)
             ENDDO
             DO im = 1,memory-1
               rho(im) = rho(im+1)
             ENDDO
             change_in_F(memory,:,:) = force-Fold
             change_in_R(memory,:,:) = R-Rold
             ! apply periodic boundary conditions to these vector differences
             CALL set_pbc(change_in_R(memory,:,:))
             rho(memory) = 1.0_q/SUM(change_in_F(memory,:,:)*change_in_R(memory,:,:))
            ENDIF
        ELSE
           ! reset down the force
           itr = 0
        ENDIF

        Rold = R
        Fold = force
        ! compute Ho*g
        IF (itr .LT. memory) THEN
          bound = itr
        ELSE
          bound = memory
        ENDIF
        direction = force
        IF (itr .NE. 0) THEN
          DO im=1,bound 
            jm = (bound+1-im)
            alpha(im) = rho(jm)*SUM(change_in_R(jm,:,:)*direction)
            direction = direction-(alpha(im)*change_in_F(jm,:,:))
          ENDDO
          ! direction = Ho*direction !Ho is idenity matrix
          DO im=1,bound
            beta = rho(im)*SUM(change_in_F(im,:,:)*direction)
            direction = direction+change_in_R(im,:,:)*(alpha(im)-beta)
          ENDDO
        ENDIF
        direction =direction/SQRT(SUM(direction*direction))

        ! finite step down force
        R = R + direction*finite_step
      ELSE
        fdstep = .TRUE.
        optflag = .FALSE. ! gives control back to method
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
            step_size=SIGN(maxmove,step_size)
          ELSE
            step_size = step_size-0.5_q*finite_step  !(*)
          ENDIF
        ENDIF
! Move now from the configuration after the fd_step, so (*) has a "-" sign
        R = R+(direction*step_size)
      ENDIF

      ! convert position back to direct cord
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

    END SUBROUTINE lbfgs_step

!**********************************************************************
!
! lbfgs init
!
!**********************************************************************

    SUBROUTINE lbfgs_init(T_INFO,IO)
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

      ! read in varables used for lbfgs

      memory=4
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'NLBFGSMEM','=','#',';','I', &
     &            memory,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''NBFGSMEM'' from file INCAR.'
         STOP
      ENDIF

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

      ! initialize the variables: vectors and matricies
      ALLOCATE(change_in_F(memory,3,nions),change_in_R(memory,3,nions))
      ALLOCATE(R(3,nions),Fold(3,nions),Rold(3,nions),direction(3,nions))
      ALLOCATE(rho(memory),alpha(memory))


      fdstep=.TRUE.
    END SUBROUTINE lbfgs_init

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
      REAL(q) :: v1(3,nions),v2(3,nions),vproj(3,nions)
      vproj=v2*SUM(v1*v2)/SQRT(SUM(v2*v2))
      END FUNCTION vproj


  END MODULE lbfgs
