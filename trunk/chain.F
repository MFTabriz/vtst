#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: chain.F,v 1.9 2006-09-28 20:27:14 graeme Exp $
!
! Module which controlls the running of four methods, the nudged
! elastic band (neb.F), the dimer (dimer.F), lanczos (lanczos.F), and
! the dynamical matrix method (dynmat.F).  The purpose of this module
! is to determine which of these methods should be run, and to do so at
! each ionic step.
!
! TODO: add a module to separate optimizers from the transition
! state methods.  The conjugate gradient and quick-min algorithms
! will also work with normal minimization runs.
!
! NOTE: the vasp folks have now added their own dynamical matrix
! code.  The only advantage of this dynamical matrix implementation
! is that you can combine forces from multiple vasp runs.  This
! makes it easier to separate a large calculation into serveral jobs,
! and systematiCALLy check for convergence of the normal modes, or
! prefactors of reactions.
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
    USE opt

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chain_force, chain_init
    PUBLIC :: Sum_Chain,And_Chain,LHYPER_NUDGE
    INTEGER ICHAIN
!    REAL SPRING
    LOGICAL optflag

!**********************************************************************
!
!  General force routine for any method using the repeated image mode.
!  The variable ICHAIN determines which method to use (see Chain_Init)
!
!**********************************************************************

  CONTAINS
    SUBROUTINE chain_force(nions,posion,toten,force,a,b,iu6)
      INTEGER :: nions,ni,nj
      INTEGER :: iu6
      REAL(q) :: ftotmag, fmaxatom,ftemp
      REAL(q) :: posion(3,nions)
      REAL(q) :: force(3,nions),toten
      REAL(q) :: a(3,3),b(3,3)
      LOGICAL :: optflag

      ! optflag indicates who has control
      !  true: optimizer is active
      !  false: chain method is active (for dimer/lanczos)
      IF(IMAGES==0) THEN
        optflag=.TRUE.
      ELSE
        IF(ICHAIN==0) THEN
          CALL neb_force(optflag,nions,posion,toten,force,a,b,iu6)
        ELSEIF(ICHAIN==1) THEN
          CALL dynmat_force(optflag,nions,posion,toten,force,a,b,iu6)
        ELSEIF(ICHAIN==2) THEN
          CALL dimer_force(optflag,nions,posion,toten,force,a,b,iu6)
        ELSEIF(ICHAIN==3) THEN
          CALL lanczos_force(optflag,nions,posion,toten,force,a,b,iu6)
        ENDIF
      ENDIF

      ftotmag=0.
      fmaxatom=0.
      DO ni=1,nions
         ftemp=0.
         DO nj=1,3
            ftotmag=ftotmag+force(nj,ni)**2
            ftemp=ftemp+force(nj,ni)**2
         ENDDO
         IF (ftemp.gt.fmaxatom) fmaxatom=ftemp
      ENDDO
      ftotmag=SQRT(ftotmag/REAL(nions))
      fmaxatom=SQRT(fmaxatom)
      IF (iu6>=0) WRITE(iu6,4693) fmaxatom,ftotmag
 4693 format(1x,' FORCES: max atom, RMS ',2f12.6)

      ! our own optimizers (optflag: do or do not optimize)  
      CALL opt_step(optflag,nions,posion,toten,force,a,b,iu6)

    END SUBROUTINE chain_force

!**********************************************************************
!
! Initialize the chain (repeated image mode) and determine which of the
! three possible methods to use based on the ICHAIN variable:
!   ICHAIN==0: nudged elastic band (default)
!   ICHAIN==1: dynamical matrix
!   ICHAIN==2: dimer method
!   ICHAIN==3: lanczos method
!
!**********************************************************************

    SUBROUTINE chain_init (T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO

      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt) :: LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics) :: DYN

      INTEGER :: IERR,N,IOPT
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM

      ! initialize optimizer
      CALL opt_init(T_INFO, IO)
      optflag=.FALSE.

      ! initialize chain based method
      ICHAIN=0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ICHAIN','=','#',';','I', &
     &            ICHAIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(IO%IU0,*)'Error reading item ''ICHAIN'' from file INCAR.'
         STOP
      ENDIF

      IOPT=0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IOPT','=','#',';','I', &
     &            IOPT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(IO%IU0,*)'Error reading item ''IOPT'' from file INCAR.'
         STOP
      ENDIF


      IF(IO%IU6>=0) WRITE(IO%IU6,*) 'CHAIN: Read ICHAIN ',ICHAIN

      IF(ICHAIN==0) THEN
        IF(IMAGES>0) THEN
          IF(IO%IU6>=0) WRITE(IO%IU6,*) 'CHAIN: Running the NEB'
        ENDIF
        CALL neb_init(T_INFO,IO)
      ELSEIF(ICHAIN==1) THEN
        IF(IO%IU6>=0) WRITE(IO%IU6,*) 'CHAIN: Running the Dynamical Matrix'
        CALL dynmat_init(T_INFO,IO)
      ELSEIF(ICHAIN==2) THEN
        IF(IO%IU6>=0) WRITE(IO%IU6,*) 'CHAIN: Running the Dimer method'
        CALL dimer_init(T_INFO,IO)
        IF(iopt==0) THEN
          WRITE(IO%IU6,*) 'CHAIN: Must set IOPT>0 to use the Dimer method'
          STOP
        ENDIF
      ELSEIF(ICHAIN==3) THEN
        IF(IO%IU6>=0) WRITE(IO%IU6,*) 'CHAIN: Running the Lanczos method'
        CALL lanczos_init(T_INFO,IO)
        IF(iopt==0) THEN
          WRITE(IO%IU6,*) 'CHAIN: Must set IOPT>0 to use the Lanczos method'
          STOP
        ENDIF
      ENDIF

    END SUBROUTINE chain_init


!**********************************************************************
!
! Returns true IF hyper nudged elastic band method is used
!
!**********************************************************************
!
      FUNCTION LHYPER_NUDGE()
      LOGICAL LHYPER_NUDGE
!      IF (images==0 .OR. spring /= 0 ) THEN
      IF (images==0 ) THEN
        LHYPER_NUDGE=.FALSE.
      ELSE
        LHYPER_NUDGE=.TRUE.
      ENDIF
      END FUNCTION LHYPER_NUDGE
!
!
!**********************************************************************
!
! MPI routines
!
!**********************************************************************
!
! Sum over elements
!
      SUBROUTINE sum_chain(value)
      REAL(q) :: value
      IF (images==0) RETURN
!      IF (spring==-1000) RETURN
      CALLMPI_C(M_sum_d(comm_chain,value,1))
      END SUBROUTINE sum_chain
!
! LOGICAL sum over elements
!
      SUBROUTINE and_chain(value)
      LOGICAL :: value
      REAL(q) :: sum_all
      IF (images==0) RETURN
!      IF (spring==-1000) RETURN
      sum_all=1
      IF (value) sum_all=0
      CALLMPI_C(M_sum_d(comm_chain,sum_all,1))
      value=.TRUE.
      IF (sum_all>=1) value=.FALSE.
      END SUBROUTINE and_chain

END MODULE chain
