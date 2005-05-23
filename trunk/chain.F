#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: chain.F,v 1.2 2005-05-23 19:51:40 graeme Exp $
!
! Module which controlls the running of three methods, the nudged
! elastic band (neb.F), the dimer (dimer.F) and the dynamical matrix
! method (dynmat.F).  The purpose of this module is to determine which
! of these methods should be run, and to do so at each ionic step.
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

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: chain_force, chain_init
    PUBLIC :: Sum_Chain,And_Chain,LHYPER_NUDGE
    INTEGER ICHAIN
    REAL SPRING

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

      IF(ICHAIN==0) THEN
        call neb_force(nions,posion,toten,force,a,b,iu6)
      ELSEIF(ICHAIN==1) THEN
        call dynmat_force(nions,posion,toten,force,a,b,iu6)
      ELSEIF(ICHAIN==2) THEN
        call dimer_force(nions,posion,toten,force,a,b,iu6)
      ELSEIF(ICHAIN==3) THEN
        call lanczos_force(nions,posion,toten,force,a,b,iu6)
      ENDIF

      ftotmag=0.
      fmaxatom=0.
      DO ni=1,nions
         ftemp=0.
         DO nj=1,3
            ftotmag=ftotmag+force(nj,ni)**2
            ftemp=ftemp+force(nj,ni)**2
         ENDDO
         if (ftemp.gt.fmaxatom) fmaxatom=ftemp
      ENDDO
      ftotmag=sqrt(ftotmag/real(nions))
      fmaxatom=sqrt(fmaxatom)
      if (iu6>=0) write(iu6,4693) fmaxatom,ftotmag
 4693 format(1x,' FORCES: max atom, RMS ',2f12.6)

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

      integer :: IERR,N
      character*1 :: CHARAC
      complex(q) :: CDUM 
      logical :: LDUM
      real(q) :: RDUM

      ICHAIN=0
      call RDATAB(.True.,'INCAR',IO%IU5,'ICHAIN','=','#',';','I', &
     &            ICHAIN,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(IO%IU0,*)'Error reading item ''ICHAIN'' from file INCAR.'
         STOP
      ENDIF

      if(IO%IU6>=0) write(IO%IU6,*) 'CHAIN: Read ICHAIN ',ICHAIN

      IF(ICHAIN==0) THEN
        IF(IMAGES>0) THEN
          if(IO%IU6>=0) write(IO%IU6,*) 'CHAIN: Running the NEB'
        ENDIF
        call neb_init(T_INFO,IO)
      ELSEIF(ICHAIN==1) THEN
        if(IO%IU6>=0) write(IO%IU6,*) 'CHAIN: Running the Dynamical Matrix'
        call dynmat_init(T_INFO,IO)
      ELSEIF(ICHAIN==2) THEN
        if(IO%IU6>=0) write(IO%IU6,*) 'CHAIN: Running the Dimer method'
        call dimer_init(T_INFO,IO)
      ELSEIF(ICHAIN==3) THEN
        if(IO%IU6>=0) write(IO%IU6,*) 'CHAIN: Running the Lanczos method'
        call lanczos_init(T_INFO,IO)
      ENDIF

    END SUBROUTINE chain_init


!**********************************************************************
!
! Returns true if hyper nudged elastic band method is used
!
!**********************************************************************
!
      function LHYPER_NUDGE()
      logical LHYPER_NUDGE
      if (images==0 .OR. spring /= 0 ) then
        LHYPER_NUDGE=.FALSE.
      else
        LHYPER_NUDGE=.TRUE.
      endif
      end function LHYPER_NUDGE
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
      subroutine sum_chain(value)
      real(q) :: value
      if (images==0) return
      if (spring==-1000) return
      CALLMPI_C(M_sum_d(comm_chain,value,1))
      end subroutine sum_chain
!
! logical sum over elements
!
      subroutine and_chain(value)
      logical :: value
      real(q) :: sum
      if (images==0) return
      if (spring==-1000) return
      sum=1
      if (value) sum=0
      CALLMPI_C(M_sum_d(comm_chain,sum,1))
      value=.TRUE.
      if (sum>=1) value=.FALSE.
      end subroutine and_chain

END MODULE chain
