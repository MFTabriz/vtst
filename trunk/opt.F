#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: opt.F,v 1.15 2006-11-10 09:44:03 dano Exp $
!
! Module which implements our optimizers
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE opt
    USE prec
    USE lattice
    USE lbfgs
    USE cg
    USE qm
!    USE bfgs
    USE sd
    USE fire
!    USE dynamic
    

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: opt_step, opt_init  !call opt_init from chain_init

    INTEGER :: nions,iu6,iopt
    REAL(q) :: dir2car(3,3),car2dir(3,3)

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

      ! if optflag is false, do nothing
      IF(.NOT. optflag) RETURN

      ! things we want to do for all optimizers?

      IF (iopt==1) THEN
        CALL lbfgs_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      ELSEIF (iopt==2) THEN
        CALL cg_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      ELSEIF (iopt==3) THEN
        CALL qm_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      ELSEIF (iopt==4) THEN
        CALL sd_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
!      ELSEIF (iopt==5) THEN
!        CALL bfgs_step(optflag,tinfo_nions,posion,toten,force, &
!                        latt_a,latt_b,tinfo_iu6)
!      ELSEIF (iopt==6) THEN
!        CALL dynamic_step(optflag,tinfo_nions,posion,toten,force, &
!                        latt_a,latt_b,tinfo_iu6)
      ELSEIF (iopt==7) THEN
        CALL fire_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      ENDIF

      ! call periodic boundary conditions
      CALL set_pbc(posion)

    END SUBROUTINE opt_step

!**********************************************************************
!
! optimizer initilizer
!
!**********************************************************************

    SUBROUTINE opt_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,ibrion
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM,potim

      nions=T_INFO%nions
      iu6=IO%IU6

      ! read iopt (0=use vasp optimizers, 1=lbfgs, 2=cg, 3=qm,
      !       4=sd, 5=bfgs, 6=md, 7=fire)

      iopt = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IOPT','=','#',';','I', &
     &            iopt,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''IOPT'' from file INCAR.'
         STOP
      ENDIF

      potim = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'POTIM','=','#',';','F', &
     &            IDUM,potim,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''POTIM'' from file INCAR.'
         STOP
      ENDIF

      ibrion = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IBRION','=','#',';','I', &
     &            ibrion,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''IBRION'' from file INCAR.'
         STOP
      ENDIF

      ! write out the optimizer chosen
      IF( (iopt.EQ.0) .AND. (ibrion.EQ.1) ) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) "OPT: USING VASP's QUASI-newton optimizer"
        RETURN
      ELSEIF( (iopt.EQ.0) .AND. (ibrion.EQ.2) ) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) "OPT: USING VASP's conjugate-gradient optimizer"
        RETURN
      ELSEIF( (iopt.EQ.0) .AND. (ibrion.EQ.3) ) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) "OPT: USING VASP's Qucik-min optimizer"
        RETURN
      ELSEIF (iopt.EQ.0) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) "OPT: USING VASP OPTIMIZERS"
        RETURN
      ELSEIF (iopt.EQ.1) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: USING LBFGS optimizer'
      ELSEIF (iopt.EQ.2) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: USING conjugate-gradient optimizer'
      ELSEIF (iopt.EQ.3) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: USING Quick-min optimizer'
      ELSEIF (iopt.EQ.4) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: USING Steepest-descent optimizer'
      ELSEIF (iopt.EQ.5) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: USING Full BFGS optimizer'
      ELSEIF (iopt.EQ.6) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: Verlet dynamics'
      ELSEIF (iopt.EQ.7) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: USING FIRE optimizer'
      ENDIF

      IF(potim/=0 .OR. ibrion/=3) THEN
        IF (IO%IU0>=0) &
             WRITE(IO%IU6,*) 'OPT: Must set IBRION=3 and POTIM=0 for IOPT>0'
        STOP
      ENDIF


      ! optimizer specific initialization
      
      IF (iopt==1) THEN
        CALL lbfgs_init(T_INFO,IO)
      ELSEIF (iopt==2) THEN
        CALL cg_init(T_INFO,IO)
      ELSEIF (iopt==3) THEN
        CALL qm_init(T_INFO,IO)
      ELSEIF (iopt==4) THEN
        CALL sd_init(T_INFO,IO)
!      ELSEIF (iopt==5) THEN
!        CALL bfgs_init(T_INFO,IO)
!      ELSEIF (iopt==6) THEN
!        CALL dynamics_init(T_INFO,IO)
      ELSEIF (iopt==7) THEN
        CALL fire_init(T_INFO,IO)
      ENDIF

    END SUBROUTINE opt_init


!**********************************************************************
!
! Vector Functions
!
!**********************************************************************

!======================================================================
! Sets a vector to have the smallest length consistent the the periodic
! boundary conditions. This should REALly be changed to use the
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
      REAL(q) :: v1(3,nions)
      REAL(q),DIMENSION(3,nions) :: return_unit
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

  END MODULE opt
