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
    SAVE
    PRIVATE 
    PUBLIC :: lbfgs_step, lbfgs_init  !CALL opt_init from chain_init

    INTEGER :: nions,iu6,memory,itr
    REAL(q),ALLOCATABLE ::R(:,:),Rold(:,:),Fold(:,:)
    REAL(q),ALLOCATABLE :: rho(:),alpha(:)
    REAL(q),ALLOCATABLE :: direction(:,:)
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:) :: change_in_G,change_in_R

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step,maxmove
    LOGICAL :: fdstep,reset_flag

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
      INTEGER :: bound,im,jm

      R = posion
      dir2car=LATT_A
      car2dir=LATT_B

      ! Convert the position into Cartesian coordinates
      CALL dirkar(nions,R,latt_a)

      IF (fdstep) THEN  ! Finite difference step

        fdstep = .FALSE.
        optflag = .TRUE. ! keeps contol so it can do finite step

        ! check for reset of direction
        a1 = ABS(SUM(force*Fold))
        a2 = SUM(Fold*Fold)
!        IF ((a1 .LE. 0.5_q*a2) .AND. (a2 .NE. 0.0_q)) THEN
        IF ((a1 .GT. 0.5_q*a2) .OR. (a2 .EQ. 0.0_q)) THEN
           reset_flag=.TRUE.
        ENDIF

        ! Reset the lbfgs memory
        IF(reset_flag) THEN
           IF(iu6>0) write(iu6,*) ' LBFGS, init'
           itr=0
           reset_flag=.FALSE.
        ELSE
           ! find new direction
           itr = itr+1
           IF (itr .LT. memory) THEN
             IF(iu6>0) write(iu6,*) ' LBFGS, adding new vector'

             change_in_G(itr,:,:) = -(force-Fold)
             IF(iu6>0) write(iu6,*) ' LBFGS, DeltaG:'
             IF(iu6>0) write(iu6,'(3f16.8)') change_in_G(itr,:,:)

             change_in_R(itr,:,:) = R-Rold 
             ! apply periodic boundary conditions to these vector differences
             CALL set_pbc(change_in_R(itr,:,:))
             IF(iu6>0) write(iu6,*) ' LBFGS, DeltaR:'
             IF(iu6>0) write(iu6,'(3f16.8)') change_in_R(itr,:,:)

             rho(itr) = 1.0_q/SUM(change_in_G(itr,:,:)*change_in_R(itr,:,:))
             IF(iu6>0) write(iu6,*) ' LBFGS, itr,rho: ',itr,rho(itr)
           ELSE
             IF(iu6>0) write(iu6,*) ' LBFGS, shifting vectors'
             DO im = 1,memory-1
               change_in_G(im,:,:) = change_in_G(im+1,:,:)
             ENDDO
             DO im = 1,memory-1
               change_in_R(im,:,:) = change_in_R(im+1,:,:)
             ENDDO
             DO im = 1,memory-1
               rho(im) = rho(im+1)
!               IF(iu6>0) write(iu6,*) ' LBFGS, itr,rho: ',itr,rho(itr)
             ENDDO
             change_in_G(memory,:,:) = -(force-Fold)
             change_in_R(memory,:,:) = R-Rold
             ! apply periodic boundary conditions to these vector differences
             CALL set_pbc(change_in_R(memory,:,:))
             rho(memory) = 1.0_q/SUM(change_in_G(memory,:,:)*change_in_R(memory,:,:))
            ENDIF
        ENDIF

        Rold = R
        Fold = force

        ! compute Ho*g
        IF (itr .LT. memory) THEN
          bound = itr
        ELSE
          bound = memory
        ENDIF
        IF(iu6>0) write(iu6,*) ' LBFGS, bound: ',bound
        direction = -1._q*force
        
        ! calculate Ho*g
        IF (itr .NE. 0) THEN
          DO im=1,bound 
            jm = (bound+1-im)
            alpha(jm) = rho(jm)*SUM(change_in_R(jm,:,:)*direction)
            direction = direction-(alpha(jm)*change_in_G(jm,:,:))
          ENDDO
          ! direction = Ho*direction !Ho is identity matrix
          DO im=1,bound
            beta = rho(im)*SUM(change_in_G(im,:,:)*direction)
            direction = direction+change_in_R(im,:,:)*(alpha(im)-beta)
          ENDDO
        ENDIF
        ! direction down gradent
        direction= -1._q*direction
        direction = direction/SQRT(SUM(direction*direction))

        IF(iu6>0) write(iu6,*) ' LBFGS, direction:'
        IF(iu6>0) write(iu6,'(3f16.8)') direction(:,:)


        ! finite step down force
        R = R + direction*finite_step

      ELSE  ! Translation step

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
          IF (ABS(step_size) .GT. maxmove) THEN
            step_size = SIGN(maxmove,step_size)-SIGN(finite_step,step_size)
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
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ILBFGSMEM','=','#',';','I', &
     &            memory,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''IBFGSMEM'' from file INCAR.'
         STOP
      ENDIF

      finite_step=0.001_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FDSTEP','=','#',';','F', &
     &            IDUM,finite_step,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''FDSTEP'' from file INCAR.'
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
      ALLOCATE(change_in_G(memory,3,nions),change_in_R(memory,3,nions))
      ALLOCATE(R(3,nions),Fold(3,nions),Rold(3,nions),direction(3,nions))
      ALLOCATE(rho(memory),alpha(memory))

      change_in_G=0._q
      change_in_R=0._q
      Fold=0._q

      reset_flag=.TRUE.
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

  END MODULE lbfgs

