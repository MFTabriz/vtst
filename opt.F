#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: opt.F,v 1.20 2007-05-07 06:44:14 graeme Exp $
!
! Module which implements forced based optimizers
!
!**********************************************************************

  MODULE opt
    USE prec
    USE lattice
    USE lbfgs
    USE cg
    USE qm
    USE bfgs
    USE sd
    USE fire
!    USE dynamic

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: opt_step, opt_init  ! call opt_init from chain_init

    INTEGER :: nions,iopt,iu6

!**********************************************************************
! Wrapper function to call the appropriate optimizer
!**********************************************************************

  CONTAINS
    SUBROUTINE opt_step(optflag,posion,toten,force,hstress,latt_a,latt_b)
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: hstress(3,3),latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      ! if optflag is false, do nothing
      IF (iu6>=0) WRITE(iu6,*) ''
      IF (iu6>=0) WRITE(iu6,*) 'OPT: Flag',optflag
      IF (.NOT. optflag) RETURN

      IF (iopt==1) THEN
        IF(iu6>=0) WRITE(iu6,*) 'OPT: LBFGS Step'
        CALL lbfgs_step(optflag,posion,toten,force,latt_a,latt_b)
      ELSEIF (iopt==2) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: CG Step'
        CALL cg_step(optflag,posion,toten,force,latt_a,latt_b)
      ELSEIF (iopt==3) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: QM Step'
        CALL qm_step(optflag,posion,toten,force,hstress,latt_a,latt_b)
      ELSEIF (iopt==4) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: SD Step'
        CALL sd_step(optflag,posion,toten,force,latt_a,latt_b)
      ELSEIF (iopt==5) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: BFGS Step'
        CALL bfgs_step(optflag,posion,toten,force,latt_a,latt_b)
!      ELSEIF (iopt==6) THEN
!        IF (iu6>=0) WRITE(iu6,*) 'OPT: Dynamic Step'
!        CALL dynamic_step(optflag,posion,toten,force,latt_a,latt_b)
      ELSEIF (iopt==7) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: FIRE Step'
        CALL fire_step(optflag,posion,toten,force,hstress,latt_a,latt_b)
      ENDIF

! Put atoms back inside box
      posion(:,:)=MOD(posion(:,:)+100._q,1._q)

    END SUBROUTINE opt_step

!**********************************************************************
! optimizer initializer
!**********************************************************************

    SUBROUTINE opt_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,ibrion,iu0
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM,potim

      iu0=IO%IU0
      iu6=IO%IU6

      ! read iopt (0=use vasp optimizers,
      !    1=lbfgs, 2=cg, 3=qm, 4=sd, 5=bfgs, 6=md, 7=fire)

      iopt = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IOPT','=','#',';','I', &
     &            iopt,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''IOPT'' from file INCAR.'
         STOP
      ENDIF

      potim = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'POTIM','=','#',';','F', &
     &            IDUM,potim,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''POTIM'' from file INCAR.'
         STOP
      ENDIF

      ibrion = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IBRION','=','#',';','I', &
     &            ibrion,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''IBRION'' from file INCAR.'
         STOP
      ENDIF

      ! write out the optimizer chosen
      IF (iu6>=0) WRITE(iu6,*) ''
      IF( (iopt.EQ.0) .AND. (ibrion.EQ.1) ) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using VASP QUASI-newton optimizer'
        RETURN
      ELSEIF( (iopt.EQ.0) .AND. (ibrion.EQ.2) ) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using VASP Conjugate-Gradient optimizer'
        RETURN
      ELSEIF( (iopt.EQ.0) .AND. (ibrion.EQ.3) ) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using VASP Quick-Min optimizer'
        RETURN
      ELSEIF (iopt.EQ.0) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using VASP Dynamics algorithm'
        RETURN
      ELSEIF (iopt.EQ.1) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using LBFGS optimizer'
      ELSEIF (iopt.EQ.2) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using Conjugate-Gradient optimizer'
      ELSEIF (iopt.EQ.3) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using Quick-Min optimizer'
      ELSEIF (iopt.EQ.4) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using Steepest-descent optimizer'
      ELSEIF (iopt.EQ.5) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using Full BFGS optimizer'
      ELSEIF (iopt.EQ.6) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Verlet Dynamics'
      ELSEIF (iopt.EQ.7) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Using FIRE optimizer'
      ENDIF

      IF(potim/=0 .OR. ibrion/=3) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: Must set IBRION=3 and POTIM=0 for IOPT>0'
        STOP
      ENDIF

      ! optimizer specific initialization
      
      IF (iopt==1) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: LBFGS, Init'
        CALL lbfgs_init(T_INFO,IO)
      ELSEIF (iopt==2) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: CG, Init'
        CALL cg_init(T_INFO,IO)
      ELSEIF (iopt==3) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: QM, Init'
        CALL qm_init(T_INFO,IO)
      ELSEIF (iopt==4) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: SD, Init'
        CALL sd_init(T_INFO,IO)
      ELSEIF (iopt==5) THEN
         IF (iu6>=0) WRITE(iu6,*) 'OPT: BFGS, Init'
        CALL bfgs_init(T_INFO,IO)
!      ELSEIF (iopt==6) THEN
!         IF (iu6>=0) WRITE(iu6,*) 'OPT: Dynamic, Init'
!        CALL dynamics_init(T_INFO,IO)
      ELSEIF (iopt==7) THEN
        IF (iu6>=0) WRITE(iu6,*) 'OPT: FIRE, Init'
        CALL fire_init(T_INFO,IO)
      ENDIF
      IF (iu6>=0) WRITE(iu6,*) ''

    END SUBROUTINE opt_init

  END MODULE opt
