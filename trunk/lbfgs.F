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
    private 
    public :: lbfgs_step, lbfgs_init  !call opt_init from chain_init

    INTEGER :: nions,iu6
    REAL(q)  R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)
    REAL(q),ALLOCATABLE :: rho(:),alpha(:)
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:) :: change_in_F,change_in_R

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step,maxmove,memory
    LOGICAL :: fdstep

!**********************************************************************
!
! Limited-memory bfgs method
!
!**********************************************************************

    SUBROUTINE lbfgs_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      REAL(q) :: beta
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,Favg,step_size
      INTEGER :: itr,bound,im,jm

      F = force
      R = posion

      IF (fdstep) THEN
        fdstep = .false.
        optflag = .true. ! keeps contol so it can do finite step
        ! check for reset of direction
        a1 = ABS(SUM(F*Fold))
        a2 = SUM(Fold*Fold)
        IF (a1 .LT. 0.5_q*a2 .AND. a2 .ne. 0.0_q) THEN
           ! find new direct
           IF (itr .LT. memory) THEN
             change_in_F(itr,:,:) = F-Fold
             change_in_R(itr,:,:) = R-Rold  ! need boundary conditions
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
             change_in_F(memory,:,:) = F-Fold
             change_in_R(memory,:,:) = R-Rold  ! need boundary conditions
             rho(memory) = 1.0_q/SUM(change_in_F(memory,:,:)*change_in_R(memory,:,:))
           itr = itr+1
        ELSE
           ! reset down the force
           itr = 1
           ! ? reset the change_in_pos/force vectors?
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

    END SUBROUTINE lbfgs_step

!**********************************************************************
!
! lbfgs init
!
!**********************************************************************

    SUBROUTINE lbfgs_init(T_INFO,IO)
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
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'NLBFGSMEM','=','#',';','F', &
     &            memory,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''NBFGSMEM'' from file INCAR.'
         STOP
      ENDIF

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

      ! initialize the variables: vectors and matricies
      ALLOCATE(change_in_F(memory,3,Nions),change_in_R(memory,3,Nions))
      ALLOCATE(F(3,Nions),R(3,Nions))
      ALLOCATE(Fold(3,Nions),Rold(3,Nions))
      ALLOCATE(rho(memory),alpha(memory))



      fdstep=.true.
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
      subroutine set_pbc(v1)
      real(q) :: v1(3,nions)
      call kardir(nions,v1,car2dir)
      v1=mod(v1+100.5_q,1._q)-0.5_q
      call dirkar(nions,v1,dir2car)
      end subroutine set_pbc
!======================================================================
! Returns a unit vector along v1
!======================================================================
      function return_unit(V1)
      real(q) :: v1(3,nions)
      real(q),dimension(3,nions) :: return_unit
      return_unit=v1*(1._q/sqrt(sum(v1*v1)))
      end function return_unit
!======================================================================
! Sets V1 to be a unit vector
!======================================================================
      subroutine set_unit(V1)
      real(q) :: v1(3,nions)
      v1=return_unit(v1)
      end subroutine set_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      function vproj(v1,v2)
      real(q) :: v1(3,nions),v2(3,nions),proj(3,nions)
      proj=v2*sum(v1*v2)/sum(v2*v2)
      end function vproj

  END MODULE lbfgs
