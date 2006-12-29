#include "symbol.inc"
#define DIMER
!**********************************************************************
! RCS:  $Id: dimer.F,v 1.18 2006-12-29 23:30:06 graeme Exp $
!
! This module implements the dimer method.  For more information refer
! to the web page:
!   http://theory.cm.utexas.edu/vtsttools/
! and the articles 
!   Henkelman and Jonsson, JCP 111, 7010 (1999)
!   Olsen, Kroes, Henkelman, Arnaldsson, and Jonsson, JCP 121, 9776 (2004).
!   Heyden, Bell, and Keil, JCP 123, 224101 (2005).
!
! Version Notes:
! 
! v2.00, Dec 2006
! - single image version of the dimer method
! - changed output file name to DIMCAR
! - added a minimum rotational force variable
! - use a 45 deg rotation (Heyden et al, JCP 123, 224101, 2005)
! - incorporate the new optimizers
!
! Graeme Henkelman
! henkelman@mail.utexas.edu
!**********************************************************************

    MODULE dimer
      USE prec
      USE main_mpi
      USE poscar
      USE lattice
      USE constant

      IMPLICIT NONE
      SAVE 
      PRIVATE
      PUBLIC :: Dimer_Force,Dimer_Init

      TYPE(in_struct) :: DIO
      TYPE(TYPE_info) :: DINFO
      INTEGER :: Nions,IU0,IU6,IUout,IUdim
      INTEGER :: i,j,MaxRot,NumRot,Seed,Itr
      REAL(q),ALLOCATABLE,DIMENSION(:,:) :: N,R,R0,F,F0,F1,F2
      REAL(q),ALLOCATABLE,DIMENSION(:,:) :: FN,FNold,GN,GNu
      REAL(q),DIMENSION(3,3) :: A,B
      REAL(q) :: U,U0,CN,Cth,dR,MaxFN,MinFN,F0r,FNr
      REAL(q) :: FN1,FN2,GamN,dTh,Th
      LOGICAL :: InitFlag,FdFlag,ModecarFlag,CGInitFlag,FdStep

    CONTAINS

!**********************************************************************
! Initialize the dimer
!**********************************************************************

    SUBROUTINE Dimer_Init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO
      EXTERNAL RMARIN

      DIO = IO
      DINFO = T_INFO
      IU6 = DIO%IU6
      IU0 = DIO%IU0
      Nions = T_info%nions
      IUout = 44
     #if defined(MPI_CHAIN)
        IF (IU6>=0) WRITE(IU0,*) 'The IMAGES tag should not be set for this implementation of the dimer method.'
        STOP
     #endif
      IF (IU6>=0) OPEN(UNIT=IUout, FILE='DIMCAR')
      IF (IU6>=0) WRITE(IUout,'(5A16)') 'F','FN','U','CN','Th'
      Itr = 0
      ALLOCATE(N(3,Nions),R(3,Nions),R0(3,Nions))
      ALLOCATE(F(3,Nions),F0(3,Nions),F1(3,Nions),F2(3,Nions))
      ALLOCATE(FN(3,Nions),FNold(3,Nions),GN(3,Nions),GNu(3,Nions))
      IF (Seed/=0) CALL RMARIN(Seed,Seed)
      CALL ReadDimerVar()
      IF (IU6>=0) CALL WriteDimerVar(IU6)
      InitFlag = .true.
    END SUBROUTINE Dimer_Init

!**********************************************************************
!  Dimer method: follow lowest curvature mode up to a saddle point
!**********************************************************************

    SUBROUTINE Dimer_Force(OptFlag,TINFO_NIONS,POSION,TOTEN,TIFOR,LATT_A,LATT_B,IO_IU6)
      LOGICAL :: OptFlag
      INTEGER :: TINFO_NIONS,IO_IU6
      REAL(q),DIMENSION(3,NIONS) :: POSION,TIFOR
      REAL(q),DIMENSION(3,3) :: LATT_A,LATT_B
      REAL(q) :: TOTEN

      U = TOTEN
      F = TIFOR
      R = POSION
      A = LATT_A
      B = LATT_B
      CALL DIRKAR(NIONS,R,A)  ! Convert position to Cartesian
      CALL SetConstraints(F)

      Itr = Itr+1
      IF (InitFlag) THEN
        IF (IU6>=0) WRITE(IU6,'(A7,A13)') 'Dimer:','Init'
        CALL InitMode()
        OptFlag = .true.
        CGInitFlag = .true.
      END IF
      IF (OptFlag) THEN  ! returned from optimizer; new position
        U0 = U
        R0 = R
        F0 = F
        R = R0+N*dR
        ! Set up for dimer rotation
        OptFlag = .false.
        FdStep = .true.
        NumRot = 0
      ELSE
        ! get rotational force
        CALL UpdateDimer()
        IF (IU6>=0) CALL WriteDimerData(IU6)
        ! check to see if rotation is converged
        IF (FNr<MinFN) THEN
          IF (IU6>=0) WRITE(IU6,'(A7,A20)') 'Dimer:','Rotation converged'
          OptFlag = .true.
        ELSE
          CALL RotateDimer()
          IF (NumRot>=MaxRot) OptFlag = .true.
        END IF
      END IF

      IF (OptFlag) THEN  ! going back to the optimizer; set R and F
        R=R0
        IF (CN<0) THEN
          F = F0-N*SUM(F0*N)*2._q
        ELSE
          F = N*SUM(F0*N)*(-1._q)
        END IF
        FdStep = .true.
        NumRot = 0
      END IF

      ! Send the current position and force back to vasp
      POSION = R
      TIFOR = F
      CALL KARDIR(NIONS,POSION,B)
      IF (IU6>=0) CALL WFORCE(IUout)  ! Empty the file buffer

    END SUBROUTINE Dimer_Force

!**********************************************************************
! Dimer Functions
!**********************************************************************

! Read or generate the initial lowest mode

    SUBROUTINE InitMode()
      REAL(q),EXTERNAL :: RANG
      ! If the MODECAR file exists, read N, otherwise use random N
      INQUIRE(FILE='MODECAR',EXIST=ModecarFlag)
      IF (ModecarFlag) THEN
        IF (IU6 > 0) WRITE(IUdim,'(A)') 'MODECAR found ... reading it'
        OPEN(210,FILE='MODECAR',ACTION='read',STATUS='old')
        READ(210,*) (N(1:3,i) , i=1,Nions)
        CLOSE(210)
      ELSE
        DO i = 1,Nions
          DO j = 1,3
            N(j,i) = RANG(0._q,1._q)
          END DO
        END DO
      END IF
      CALL SetConstraints(N)
      CALL SetUnit(N)
    END SUBROUTINE InitMode

! Update the forces and energy of the dimer

    SUBROUTINE UpdateDimer()
      F1 = F
      F2 = 2._q*F0-F1
      FN = ((F1-N*SUM(F1*N))-(F2-N*SUM(F2*N))/(2._q*dR)
      FNr = SQRT(SUM(FN**2))
      F0r = SQRT(SUM(F0**2))
      CN = (SUM(F2*N)-SUM(F1*N))/(2._q*dR)
    END SUBROUTINE UpdateDimer

! Dimer rotation

    SUBROUTINE RotateDimer()
      REAL(q) a1,a2
      ! Conjugate gradient and modified Newtons method for rotation
      IF (FdStep) THEN
        ! Finite difference step
        IF (CGInitFlag) THEN
          CGInitFlag=.false.
          FNold = 0.0_q
          GN = 0.0_q
        END IF
        a1 = ABS(SUM(FN*FNold))
        a2 = SUM(FNold*FNold)
        IF ((a1<=0.5_q*a2) .AND. (a2/=0._q)) THEN
          GamN = SUM(FN*(FN-FNold))/a2
        ELSE
          GamN = 0.0_q
        ENDIF
        GN = FN+GN*GamN
        GNu = GN/SQRT(SUM(GN*GN))
        FN1 = SUM(FN*GNu)
        CALL Rotate(N,GNu,PI/4._q)
        R = R0+N*dR
        FdStep = .false.
        IF (IU6>=0) WRITE(IU6,'(A7,A20,F14.6)') 'Dimer:','FN1',FN1
      ELSE
        ! Rotation step
        NumRot = NumRot+1
        FN2 = SUM(FN*GNu)
        IF (FN2/=0._q) THEN
          Th = ATAN(FN1/FN2)/(-2._q) - PI/4._q
        ELSE
          Th = PI/(-2._q) - PI/4._q
        END IF
        IF (FN2>0._q) Th = Th + PI/2._q
        CALL Rotate(N,GNu,Th)
        CALL SetUnit(N)
        R = R0+N*dR
        FdStep = .true.
        IF (IU6>=0) THEN
          WRITE(IU6,'(A7,A10,F14.6)') 'Dimer:','FN2',FN2
          WRITE(IU6,'(A7,A10,F14.6,/)') 'Dimer:','Th',Th
          WRITE(IUout,'(5F16.7)') F0r,FNr,U0,CN,Th
        END IF
      END IF
    END SUBROUTINE RotateDimer

!**********************************************************************
! Vector Functions
!**********************************************************************

! Sets a vector to have the smallest length consistent the the periodic boundaries
! This should really be changed to use the Wigner-Sitz cell

    SUBROUTINE SetPBC(V)
      REAL(q),DIMENSION(3,Nions) :: V
      CALL KARDIR(Nions,V,B)
      V = MOD(V+100.5_q,1._q)-0.5_q
      CALL DIRKAR(Nions,V,A)
    END SUBROUTINE SetPBC

! Rotates both V1 and V2

    SUBROUTINE Rotate(V1,V2,tTh)
      REAL(q),DIMENSION(:,:) :: V1,V2,V1tmp
      REAL(q) :: tTh,cTh,sTh
      cTh = COS(tTh)
      sTh = SIN(tTh)
      V1tmp = V1
      V1 = V1*cth+V2*sTh
      V2 = V2*cth-V1tmp*sTh
    END SUBROUTINE Rotate

! Sets V to be a unit vector

    SUBROUTINE SetUnit(V)
      REAL(q),DIMENSION(:,:) :: V
      V = V*(1._q/SQRT(SUM(V*V)))
    END SUBROUTINE SetUnit

! Set the constained coordinates in a vector to zero

    SUBROUTINE SetConstraints(V)
      ! check if the constrains are cartesian or direct
      ! also, remove drift if no constraints
      REAL(q),DIMENSION(3,Nions) :: V
      REAL(q) Vtmp(3)
      IF (DINFO%LSDYN) THEN
        DO i = 1,Nions
          Vtmp(1:3) = V(1:3,i)
          CALL KARDIR(1,Vtmp,B)
          DO j = 1,3
            IF (.NOT. DINFO%LSFOR(j,i)) Vtmp(j) = 0._q
          END DO
          CALL DIRKAR(1,Vtmp,A)
          V(1:3,i) = Vtmp(1:3)
        END DO
      ELSE
        Vtmp = SUM(V,DIM=2)
        DO i = 1,Nions
          V(1:3,i) = V(1:3,i)-Vtmp(1:3)
        END DO
        IF (IU6>=0) WRITE(IUdim,'(A20,3F14.6)') 'Dimer: Drift',Vtmp(1:3)
      END IF
    END SUBROUTINE SetConstraints

!**********************************************************************
! IO routines
!**********************************************************************

! Write initial dimer variables

    SUBROUTINE WriteDimerVar(IU)
      INTEGER IU
      WRITE(IU,'(/,A23,/)') 'Dimer: Input Parameters'
      WRITE(IU,'(A7,A20,I7)') 'Dimer:','Seed',Seed
      WRITE(IU,'(A7,A20,I7)') 'Dimer:','MaxRot',MaxRot
      WRITE(IU,'(A7,A20,F14.6)') 'Dimer:','dR',dR
      WRITE(IU,'(A7,A20,F14.6)') 'Dimer:','MaxFN',MaxFN
      WRITE(IU,'(A7,A20,F14.6,/)') 'Dimer:','MinFN',MinFN
    END SUBROUTINE WriteDimerVar

! Write dimer parameters

    SUBROUTINE WriteDimerData(IU)
      INTEGER IU
      WRITE(IU,'(A7,A13,F14.6)') 'Dimer:','CN',CN
      WRITE(IU,'(A7,A13,F14.6)') 'Dimer:','F0',F0r
      WRITE(IU,'(A7,A13,F14.6,/)') 'Dimer:','FN',FNr
    END SUBROUTINE WriteDimerData

! Read Dimer Variables from the INCAR file

    SUBROUTINE ReadDimerVar()
      INTEGER :: IDUM,IERR,INint
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM
      IF (IU6>=0) WRITE(IU6,'(A7,A20)') 'Dimer:','Loading variables'

      Seed=0
      CALL RDATAB(.true.,'INCAR',DIO%IU5,'DSeed','=','#',';','I', &
     &            Seed,RDUM,CDUM,LDUM,CHARAC,INint,1,IERR)
      MaxRot=1
      CALL RDATAB(.true.,'INCAR',DIO%IU5,'DMaxRot','=','#',';','I', &
     &            MaxRot,RDUM,CDUM,LDUM,CHARAC,INint,1,IERR)
      dR=5E-3_q
      CALL RDATAB(.true.,'INCAR',DIO%IU5,'DdR','=','#',';','F', &
     &            IDUM,dR,CDUM,LDUM,CHARAC,INint,1,IERR)
      MaxFN=1._q
      CALL RDATAB(.true.,'INCAR',DIO%IU5,'DMaxFN','=','#',';','I', &
     &            IDUM,MaxFN,CDUM,LDUM,CHARAC,INint,1,IERR)
      MinFN=0.01_q
      CALL RDATAB(.true.,'INCAR',DIO%IU5,'DMinFN','=','#',';','I', &
     &            IDUM,MinFN,CDUM,LDUM,CHARAC,INint,1,IERR)
    END SUBROUTINE ReadDimerVar

    END MODULE dimer

