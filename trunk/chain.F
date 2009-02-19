#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: chain.F,v 1.39 2009-02-19 03:44:26 graeme Exp $
!
! Module which controls the running of four methods, the nudged
! elastic band (neb.F), the dimer (dimer.F), lanczos (lanczos.F), and
! the dynamical matrix method (dynmat.F).  The purpose of this module
! is to determine which of these methods should be run, and to do so at
! each ionic step.
!
! A set of force-based optimers are also included for use with normal
! optimizations, or transition state calculations.  The optimizers are
! steepest-descent, quick-min, conjugate-gradients, and LBFGS.
!
! NOTE: the vasp folks have now added their own dynamical matrix
! code.  The only advantage of this dynamical matrix implementation
! is that you can combine forces from multiple vasp runs.  This
! makes it easier to separate a large calculation into serveral jobs,
! and systematically check for convergence of the normal modes, or
! prefactors of reactions.
!
! For more information see: http://theory.cm.utexas.edu/vtsttools/
!
! Contributers:
!   Andri Arnaldsson
!   Graeme Henkelman
!   Hannes Jonsson
!   Daniel Sheppard
!   Blas Uberuaga
!   Lijun Xu
!
! Email: henkelman@mail.utexas.edu
!
! Version 2.03d, Feb 18, 2009
!
!**********************************************************************

  MODULE chain
    USE prec
    USE main_mpi
    USE poscar
    USE lattice
    USE neb
    USE dynmat
    USE dimer
    USE lanczos
!INS_BEGIN
    USE instanton
!INS_END
    USE opt

    IMPLICIT NONE
    SAVE
    PRIVATE
    PUBLIC :: chain_force, chain_init
    PUBLIC :: Sum_Chain,And_Chain,LHYPER_NUDGE
#ifdef EAM
    PRIVATE :: EAMForce
#endif
    INTEGER :: ICHAIN,IOPT
    LOGICAL :: optflag,fconverge,ftot_flag
    REAL(q) :: EDIFFG_local,ftot_val
    REAL(q),ALLOCATABLE :: Free(:,:)

#ifdef EAM
    ! Variables used for calling the EAM potential (testing)
    LOGICAL :: eaminit
    REAL(q),ALLOCATABLE :: Ream(:,:),Feam(:,:)
    REAL(q),ALLOCATABLE :: Rvasp(:,:),Fvasp(:,:)
!    REAL(q),ALLOCATABLE :: Free(:,:)
    REAL(q) :: Uvasp,Ueam
#endif

!**********************************************************************
!  General force routine for any method using the repeated image mode.
!  The variable ICHAIN determines which method to use (see chain_init)
!**********************************************************************

  CONTAINS
    SUBROUTINE chain_force(nions,posion,toten,force,a,b,iu6)
      INTEGER :: nions,ni,nj,iu6
      REAL(q) :: ftot,frms,fmaxatom,ftemp,fmaxdim,toten
      REAL(q),DIMENSION(3,nions) :: posion,force
      REAL(q),DIMENSION(3,nions) :: posion_vasp,force_vasp,force_dimlan
      REAL(q),DIMENSION(3,3) :: a,b

      ! for dimer/lanczos, save the force and posion for the vasp stop criteria
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force_vasp=force
        posion_vasp=posion
      ENDIF

#ifdef EAM
      IF (.NOT. eaminit) THEN
        eaminit=.TRUE.
        Ream=posion
        Rvasp=posion
        Fvasp=force
        Uvasp=toten
      ENDIF
      posion=Ream
      IF (iu6>0) WRITE(iu6,*) 'Calling EAMForce'
      IF (iu6>0) WRITE(iu6,*) 'posion'
      IF (iu6>0) WRITE(iu6,'(3F14.6)') posion
      CALL EAMForce(iu6,nions,posion,a,Feam,Ueam)
      force=Feam
      toten=Ueam
      IF (iu6>0) WRITE(iu6,*) 'force'
      IF (iu6>0) WRITE(iu6,'(3F14.6)') force
      IF (iu6>0) WRITE(iu6,*) 'U ',Ueam
#endif

      ! optflag indicates who has control
      !  true: optimizer is active
      !  false: chain method is active (for dimer/lanczos)
      IF (IMAGES==0 .AND. ICHAIN==0) THEN
        optflag=.TRUE.
      ELSE
        IF (ICHAIN==0) THEN
          CALL neb_step(optflag,posion,toten,force,a,b)
        ELSEIF (ICHAIN==1) THEN
          CALL dynmat_step(optflag,posion,toten,force,a,b)
        ELSEIF (ICHAIN==2) THEN
          CALL dimer_step(optflag,posion,toten,force,a,b)
          force_dimlan=force
        ELSEIF (ICHAIN==3) THEN
          CALL lanczos_step(optflag,posion,toten,force,a,b)
          force_dimlan=force
!INS_BEGIN
        ELSEIF (ICHAIN==4) THEN
          CALL instanton_step(optflag,posion,toten,force,a,b)
!INS_END
        ENDIF
      ENDIF

      ! zero out any added forces on frozen atoms
      force = force*Free

      ! for dimer/lanczos, check force criteria using the vasp force
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force_dimlan=force
        force=force_vasp
      ENDIF

      ftot=0._q
      fmaxatom=0._q
      fmaxdim=0._q
      DO ni=1,nions
        ftemp=0._q
        DO nj=1,3
          ftot=ftot+force(nj,ni)**2
          ftemp=ftemp+force(nj,ni)**2
          IF(fmaxdim.LT.ABS(force(nj,ni))) THEN
            fmaxdim=ABS(force(nj,ni))
          ENDIF
        ENDDO
        IF (ftemp .GT. fmaxatom) fmaxatom=ftemp
      ENDDO
      frms=SQRT(ftot/REAL(nions))
      fmaxatom=SQRT(fmaxatom)
      IF (iu6>=0) WRITE(iu6,4693) fmaxatom,frms
 4693 format(1x,' FORCES: max atom, RMS ',2f12.6)
      IF (iu6>=0) WRITE(iu6,4694) SQRT(ftot),fmaxdim
 4694 format(1x,' FORCE total and by dimension',2f12.6)

      ! for dimer/lanczos, use projected force in optimizer
      IF (ICHAIN==2 .OR. ICHAIN==3) THEN
        force=force_dimlan
      ENDIF

      ! stops based on the Magnitude of the Force
      IF(ftot_flag) THEN
        fconverge = (SQRT(ftot) .LT. ftot_val)
        CALLMPI_C(and_chain(fconverge))
        IF (fconverge) THEN
          IF (iu6>=0) WRITE(iu6,*) 'CONVERGED based on Magnitude of Force'
          STOP
        ENDIF
      ENDIF

      IF (IOPT .NE. 0) THEN
        fconverge = (fmaxatom .LT. ABS(EDIFFG_local))
        CALLMPI_C(and_chain(fconverge))
        IF (.NOT. fconverge) THEN
          ! our own optimizers (optflag: do or do not optimize)  
          CALL opt_step(optflag,posion,toten,force,a,b)
        ELSE
          IF (iu6>=0) WRITE(iu6,*) 'OPT: skip step - force has converged'
          ! GH: remove any finite difference steps taken by dimer or lan
          IF (ICHAIN==2 .OR. ICHAIN==3)  posion=posion_vasp
          IF (ICHAIN==2) THEN
            CALL dimer_fin()
          ENDIF
        ENDIF
      ENDIF

#ifdef EAM
      Ream=posion
      ! copy saved variables back for vasp
      posion=Rvasp
      force=Fvasp
      toten=Uvasp
#endif

    ! for dimer/lanczos, return the true force to vasp
    IF (ICHAIN==2 .OR. ICHAIN==3) THEN
      force=force_vasp
    ENDIF

    END SUBROUTINE chain_force

!**********************************************************************
! Initialize the chain (repeated image mode) and determine which of the
! three possible methods to use based on the ICHAIN variable:
!   ICHAIN==0: nudged elastic band (default)
!   ICHAIN==1: dynamical matrix
!   ICHAIN==2: dimer method
!   ICHAIN==3: lanczos method
!**********************************************************************

    SUBROUTINE chain_init(T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO
      INTEGER :: NI,NJ,IU0,IU6

      INTEGER :: IERR,N,IDUM
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM

#if defined(EAM)
      external EAM_POTINIT
#endif
      IU0=IO%IU0
      IU6=IO%IU6

      ! write the version number
      IF(IU6>=0) WRITE(IU6,'(/,A,/)') ' VTST: version 2.03d, (02/18/09)'

      IF(IU6>=0) WRITE(IU6,*) 'CHAIN: initializing optimizer'

      ! initialize optimizer
      CALL opt_init(T_INFO, IO)
      optflag=.FALSE.

      ! initialize chain based method
      ICHAIN=0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ICHAIN','=','#',';','I', &
     &            ICHAIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ICHAIN'' from file INCAR.'
        STOP
      ENDIF

      IOPT=0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IOPT','=','#',';','I', &
     &            IOPT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''IOPT'' from file INCAR.'
        STOP
      ENDIF

      EDIFFG_local=0.1_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'EDIFFG','=','#',';','F', &
     &            IDUM,EDIFFG_local,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''EDIFFG'' from file INCAR.'
        STOP
      ENDIF

!     USED to converge based on Magnitude of the force
      ftot_flag=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FMAGFLAG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ftot_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''FMAGFLAG'' from file INCAR.'
        STOP
      ENDIF

      ftot_val=0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FMAGVAL','=','#',';','F', &
     &            IDUM,ftot_val,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''FMAGVAL'' from file INCAR.'
        STOP
      ENDIF

      ! make sure that convergence is force based when using IOPT
      IF((IOPT .NE. 0) .AND. (EDIFFG_local .GT. 0.0_q)) THEN
        IF(IU6>=0) WRITE(IU6,*) 'Must set  EDIFFG < 0 when using IOPT > 0'
        STOP
      ENDIF

      IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Read ICHAIN ',ICHAIN

      IF(ICHAIN==0) THEN
        IF(IMAGES>0) THEN
          IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the NEB'
        ENDIF
        CALL neb_init(T_INFO,IO)
      ELSEIF(ICHAIN==1) THEN
        IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Dynamical Matrix'
        CALL dynmat_init(T_INFO,IO)
      ELSEIF(ICHAIN==2) THEN
        IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Dimer method'
        CALL dimer_init(T_INFO,IO)
        IF(iopt==0) THEN
          IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Must set IOPT>0 to use the Dimer method'
          STOP
        ENDIF
      ELSEIF(ICHAIN==3) THEN
        IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Lanczos method'
        CALL lanczos_init(T_INFO,IO)
        IF(iopt==0) THEN
          IF(IU6>=0) WRITE(IU6,*) 'CHAIN: Must set IOPT>0 to use the Lanczos method'
          STOP
        ENDIF
!INS_BEGIN
      ELSEIF (ICHAIN==4) THEN
        IF (IO%IU6>=0) WRITE(IU6,*) 'CHAIN: Running the Instanton method'
        CALL instanton_init(T_INFO,IO)
!INS_END
      ENDIF

      ! Make vector to zero out force on frozen atoms 
      ALLOCATE(Free(3,T_INFO%nions))
      Free(:,:)=1._q
      IF(T_INFO%LSDYN) THEN
        DO NI=1,T_INFO%nions
          DO NJ=1,3
            IF (.NOT.T_INFO%LSFOR(NJ,NI)) Free(NJ,NI)=0._q
          ENDDO
        ENDDO
      ENDIF

#if defined(EAM)
      CALL EAM_POTINIT()
      eaminit=.FALSE.
      ALLOCATE(Feam(3,T_INFO%nions),Ream(3,T_INFO%nions))
      ALLOCATE(Fvasp(3,T_INFO%nions),Rvasp(3,T_INFO%nions))
!      ALLOCATE(Free(3,T_INFO%nions))
!
!      Free(:,:)=1._q
!      IF(T_INFO%LSDYN) THEN
!        DO NI=1,T_INFO%nions
!          DO NJ=1,3
!            IF (.NOT.T_INFO%LSFOR(NJ,NI)) Free(NJ,NI)=0._q
!          ENDDO
!        ENDDO
!      ENDIF
#endif

    END SUBROUTINE chain_init

!**********************************************************************
! Returns true if hyper nudged elastic band method is used
!**********************************************************************

      FUNCTION LHYPER_NUDGE()
      LOGICAL LHYPER_NUDGE
!      IF (images==0 .OR. spring /= 0 ) THEN
      IF (images==0 ) THEN
        LHYPER_NUDGE=.FALSE.
      ELSE
        LHYPER_NUDGE=.TRUE.
      ENDIF
      END FUNCTION LHYPER_NUDGE

!**********************************************************************
! MPI routines
!**********************************************************************

! Sum over elements

      SUBROUTINE sum_chain(value)
      REAL(q) :: value
      IF (images==0) RETURN
      CALLMPI_C(M_sum_d(comm_chain,value,1))
      END SUBROUTINE sum_chain

! Logical sum over elements

      SUBROUTINE and_chain(value)
      LOGICAL :: value
      REAL(q) :: sum_all
      IF (images==0) RETURN
      sum_all=1
      IF (value) sum_all=0
      CALLMPI_C(M_sum_d(comm_chain,sum_all,1))
      value=.TRUE.
      IF (sum_all>=1) value=.FALSE.
      END SUBROUTINE and_chain

#if defined(EAM)

!**********************************************************************
! Empirical potential, for testing purposes only
!**********************************************************************

      SUBROUTINE EAMForce(iu6,NIONS,POSION,A,Feam,Ueam)
      INTEGER :: NIONS,iu6
      REAL(q) :: POSION(3,NIONS)
      REAL(q) :: Reamcar(3,NIONS),Feam(3,NIONS)
      REAL(q) :: A(3,3)
      INTEGER :: NI,NJ

      REAL*8 :: Rtmp(3*NIONS),Ftmp(3*NIONS)
      REAL*8 :: Ueam
      REAL*8 :: Ax,Ay,Az

      external EAM_FORCE

      Ax=A(1,1)
      Ay=A(2,2)
      Az=A(3,3)

      IF (IU6>0) WRITE(IU6,*) 'IN EAMForce'
      ! Convert positions to Cartesian
      Reamcar=POSION
      call DIRKAR(NIONS,Reamcar,A)
      IF (IU6>0) WRITE(IU6,*) 'Reamcar'
      IF (IU6>0) WRITE(IU6,'(3F14.8)') Reamcar(:,:)
      ! 
      do NI=1,NIONS
         do NJ=1,3
            Rtmp(3*(NI-1)+NJ)=Reamcar(NJ,NI)
         enddo
      enddo
      call EAM_FORCE(Nions,Rtmp,Ftmp,Ueam,Ax,Ay,Az)
      do NI=1,NIONS
         do NJ=1,3
            Feam(NJ,NI)=Ftmp(3*(NI-1)+NJ)
         enddo
      enddo
      Feam=Feam*Free
      IF (IU6>0) WRITE(IU6,*) 'EAMForce'
      IF (IU6>0) WRITE(IU6,'(3F14.8)') FEAM(:,:)
      IF (IU6>0) WRITE(IU6,*) 'Ueam: ',Ueam
      end subroutine EAMForce
#endif

END MODULE chain
