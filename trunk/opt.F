#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: opt.F,v 1.8 2006-06-28 19:07:57 dano Exp $
!
! Module which implements our optimizers
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE opt
    USE prec
    USE lattice

    IMPLICIT NONE
    private :: sd_step, cg_step, qm_step, lbfgs_step, bfgs_step
    private :: sd_init, cg_init, qm_init, lbfgs_init, bfgs_init

    public :: opt_step, opt_init  !call opt_init from chain_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step ! ...
    LOGICAL :: fdstep, optflag_local

!**********************************************************************
!
! Routine to optimize
!
!**********************************************************************

  CONTAINS
    SUBROUTINE opt_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      ! local variables

      nions=tinfo_nions
      iu6=tinfo_iu6
      dir2car=latt_a
      car2dir=latt_b
      optflag_local=optflag

      ! if optflag is false, do nothing
      if(.not. optflag) return

      ! things we want to do for all optimizers?

     if (iopt==1) then
        call lbfgs_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      elseif (iopt==2) then
        call cg_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      elseif (iopt==3) then
        call qm_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      elseif (iopt==4) then
        call sd_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      elseif (iopt==5) then
        call bfgs_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      endif

      ! call periodic boundary conditions
      optflag=optflag_local

    END SUBROUTINE opt_step


    SUBROUTINE opt_init(T_INFO,IO)
      type(in_struct) :: IO
      type(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,ibrion,iopt
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM,potim

      nions=T_INFO%nions
      iu6=IO%IU6

      ! read iopt (0=use vasp optimizers, 1=lbfgs, 2=cg, 3=qm, 4=bfgs, 5=sd)

      iopt=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'IOPT','=','#',';','F', &
     &            IOPT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''IOPT'' from file INCAR.'
         STOP
      ENDIF
      ! write out the optimizer chosen

      if (iopt==0) return

      potim=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'potim','=','#',';','F', &
     &            IDUM,POTIM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''POTIM'' from file INCAR.'
         STOP
      ENDIF

      ibrion=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'IBRION','=','#',';','F', &
     &            IBRION,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''IBRION'' from file INCAR.'
         STOP
      ENDIF

      IF(potim/=0 .OR. ibrion/=3) THEN
        WRITE(IO%IU6,*) 'OPT: Must set IBRION=3 and POTIM=0 for IOPT>0'
        STOP
      ENDIF

      ! allocate memory for local F and R

      allocate(F(3,Nions),R(3,Nions))
      allocate(Fold(3,Nions),Rold(3,Nions))

      ! optimizer specific initialization
      
      if (iopt==1) then
        call lbfgs_init(T_INFO,IO)
      elseif (iopt==2) then
        call cg_init(T_INFO,IO)
      elseif (iopt==3) then
        call qm_init(T_INFO,IO)
      elseif (iopt==4) then
        call sd_init(T_INFO,IO)
      elseif (iopt==5) then
        call bfgs_init(T_INFO,IO)
      endif

    end subroutine opt_init


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

END MODULE opt
