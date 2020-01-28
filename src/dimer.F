#include "symbol.inc"
!**********************************************************************
!
! This module implements the dimer method.  For more information refer
! to the web page:
!   http://theory.cm.utexas.edu/vtsttools/
! and the articles 
!   Henkelman and Jonsson, JCP 111, 7010 (1999)
!   Olsen, Kroes, Henkelman, Arnaldsson, and Jonsson, JCP 121, 9776 (2004).
!   Heyden, Bell, and Keil, JCP 123, 224101 (2005).
!   Kastner and Sherwood, JCP 128, 014106 (2008).
!
! Graeme Henkelman
! henkelman@utexas.edu
!
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
      PUBLIC :: Dimer_Step, Dimer_Init, Dimer_Fin

      TYPE(type_info) :: DINFO
      INTEGER :: Nions, IU0, IU5, IU6
      INTEGER :: i, j, RotMax, RotNum, Seed, Itr
      INTEGER,PARAMETER :: IUdim=44, IUcent=45
      REAL(q),ALLOCATABLE,DIMENSION(:,:) :: N, R, R0, F, F0, F1, F2
      REAL(q),ALLOCATABLE,DIMENSION(:,:) :: Feff, FN, FNold, GN, GNu
      REAL(q),DIMENSION(3,3) :: A, B
      REAL(q) :: U, U0, CN, Cth, dR, FNMax, FNMin, F0r, FNr
      REAL(q) :: FN1, FN2, GamN, dTh, Th, CN1, FN1r
      LOGICAL :: InitFlag, FdFlag, ModecarFlag, CGInitFlag, FdStep, NewFlag

    CONTAINS

!**********************************************************************
! Initialize the dimer
!**********************************************************************

    SUBROUTINE Dimer_Init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      EXTERNAL RMARIN

      DINFO = T_INFO
      IU0 = IO%IU0
      IU5 = IO%IU5
      IU6 = IO%IU6
      Nions = T_info%nions
#if defined(MPI_CHAIN)
        IF (IU6>=0) WRITE(IU0,*) 'The IMAGES tag should not be set for this implementation of the dimer method.'
        STOP
#endif
      IF (IU6>=0) OPEN(UNIT=IUdim, FILE='DIMCAR')
      IF (IU6>=0) WRITE(IUdim,'(A5,5A14)') 'Step', 'Force', 'Torque', 'Energy', 'Curvature', 'Angle'
      Itr = 1
      ALLOCATE(N(3,Nions),R(3,Nions),R0(3,Nions),Feff(3,Nions))
      ALLOCATE(F(3,Nions),F0(3,Nions),F1(3,Nions),F2(3,Nions))
      ALLOCATE(FN(3,Nions),FNold(3,Nions),GN(3,Nions),GNu(3,Nions))
      IF (Seed/=0) CALL RMARIN(Seed,Seed)
      CALL ReadDimerVar()
      InitFlag = .TRUE.

    END SUBROUTINE Dimer_Init

!**********************************************************************
!  Dimer method: follow lowest curvature mode up to a saddle point
!**********************************************************************

    SUBROUTINE Dimer_Step(OptFlag, POSION, TOTEN, TIFOR, LATT_A, LATT_B)
      LOGICAL :: OptFlag
      INTEGER :: TINFO_NIONS, IO_IU6
      REAL(q),DIMENSION(3,NIONS) :: POSION, TIFOR
      REAL(q),DIMENSION(3,3) :: LATT_A, LATT_B
      REAL(q) :: TOTEN

      ! Optimizer has control; find effective force and quit
      IF (OptFlag) THEN
        F0 = TIFOR
        CALL ProjectDimer()
        TIFOR = Feff
        RETURN
      END IF

      ! Dimer method has control
      U = TOTEN
      F = TIFOR
      R = POSION
      A = LATT_A
      B = LATT_B

      CALL DIRKAR(Nions, R, A)  ! Convert position to Cartesian
      CALL SetConstraints(F)

      IF (InitFlag) THEN
        CALL InitMode()
        InitFlag = .FALSE.
        NewFlag = .TRUE.
        FdStep = .TRUE.
        CGInitFlag = .TRUE.
        CN = 0._q
      END IF

      IF (IU6>=0.AND.NewFlag) WRITE(IU6,*) 'Dimer: -----------------' 
      IF (IU6>=0) WRITE(IU6,'(A,I5,A,I5,/)') ' Dimer: Itr', Itr, ' Rot', RotNum

      ! Find effective force except after a finite difference step

      ! returned from optimizer; new position
      IF (NewFlag) THEN
        U0 = U
        R0 = R
        F0 = F
        F0r = SQRT(SUM(F0**2))
        IF (IU6>=0) WRITE(IU6,'(A)') ' Dimer: Central Point'
        IF (IU6>=0) WRITE(IU6,'(A,F14.6,/)') ' Dimer: F0', F0r
        ! Write CENTCAR
        POSION = R0
        CALL KARDIR(NIONS, POSION, B)
        IF (IU6>=0) THEN
          OPEN(UNIT=IUcent, FILE='CENTCAR', STATUS='UNKNOWN')
          CALL OUTCENT(IUcent, .TRUE., DINFO%SZNAM2, DINFO, 1._q, A, DINFO%NTYP, & 
                    & DINFO%NITYP, DINFO%LSDYN, DINFO%NIONS, POSION, DINFO%LSFOR)
          CLOSE(IUcent)
        END IF
        CALL ProjectDimer()
        R = R0 + N*dR
        ! Set up for dimer rotation
        NewFlag = .FALSE.
        FdStep = .TRUE.
        RotNum = 1
      ELSE
        ! get rotational force
        CALL UpdateDimer()
        IF (FdStep) THEN
          CN = (SUM(F2*N) - SUM(F1*N))/(2._q*dR)
          CALL ProjectDimer()
          IF (IU6>=0) THEN
            WRITE(IU6,'(A,F14.6)') ' Dimer: CN ', CN
            WRITE(IU6,'(A,F14.6)') ' Dimer: F0 ', F0r
            WRITE(IU6,'(A,F14.6,/)') ' Dimer: FN ', FNr
          END IF
          !GH: moved this here so that the NEWMODECAR reflects the last N
          CALL WriteMode()
        END IF

        ! check to see if rotation is converged
        IF (FNr<FNMin) THEN
          IF (IU6>=0) WRITE(IU6,*) 'Dimer: Rotation converged'
          IF (IU6>=0) WRITE(IUdim,'(I5,5F14.5)') Itr, F0r, FN1r, U0, CN,0._q
          OptFlag = .TRUE.
        ELSE
          CALL RotateDimer()
          IF (FdStep) CALL ProjectDimer()
          IF (RotNum>RotMax) THEN
            IF(IU6>=0) WRITE(IU6,'(A,/)') ' Dimer: RotNum > RotMax'
            OptFlag = .TRUE.
          ELSE IF (FdStep.AND.FN1r<FNMax) THEN
            IF (IU6>=0) WRITE(IU6,'(A,/)') ' Dimer: FN < FNMax'
            OptFlag = .TRUE.
          END IF
        END IF
      END IF

      IF (OptFlag) THEN  ! going back to the optimizer; set R and Feff
        R=R0
        NewFlag = .TRUE.
        FdStep = .TRUE.
        RotNum = 0
        Itr = Itr + 1
      END IF

      ! Send the current position and force back to vasp
      POSION = R
      TIFOR = Feff
      CALL KARDIR(NIONS,POSION,B)
      IF (IU6>=0) CALL WFORCE(IUdim)  ! Empty the file buffer

    END SUBROUTINE Dimer_Step

!**********************************************************************
! Dimer cleanup
!**********************************************************************

! Dimer cleanup routine

    SUBROUTINE Dimer_Fin()
      IF (IU6>=0) WRITE(IUdim,'(I5,F14.5,7x,A3,4x,F14.5,7x,A3,11x,A3,4x)') &
      &  Itr,F0r,'---',U0,'---','---'
    END SUBROUTINE Dimer_Fin

!**********************************************************************
! Dimer Functions
!**********************************************************************

! Read or generate the initial lowest mode

    SUBROUTINE InitMode()
      REAL(q),EXTERNAL :: RANG
      IF (IU6>=0) WRITE(IU6,*) 'Dimer: Init Mode:'
      ! If the MODECAR file exists, read N, otherwise use random N
      INQUIRE(FILE='MODECAR', EXIST=ModecarFlag)
      IF (ModecarFlag) THEN
        IF (IU6>=0) WRITE(IU6,*) 'Dimer: MODECAR found'
        OPEN(210, FILE='MODECAR', ACTION='read', STATUS='old')
        READ(210,*) (N(1:3, i), i = 1, Nions)
        CLOSE(210)
      ELSE
        IF (IU6>=0) WRITE(IU6,*) 'Dimer: No MODECAR found'
        DO i = 1, Nions
          DO j = 1, 3
            N(j, i) = RANG(0._q, 1._q)
          END DO
        END DO
      END IF
      IF (IU6>=0) WRITE(IU6,*) 'Dimer: initial mode' 
      IF (IU6>=0) WRITE(IU6,'(3ES20.10)') (N(:,i) , i=1,Nions)
      CALL SetConstraints(N)
      CALL SetUnit(N)
      CALL WriteMode()
    END SUBROUTINE InitMode

! Update the forces and energy of the dimer

    SUBROUTINE UpdateDimer()
      F1 = F
      F2 = 2._q*F0 - F1
      FN = ((F1 - N*SUM(F1*N)) - (F2 - N*SUM(F2*N)))/(2._q*dR)
      FNr = SQRT(SUM(FN**2))
    END SUBROUTINE UpdateDimer

! Force projection to find the effective force

    SUBROUTINE ProjectDimer()
      IF (iu6>=0) THEN
        WRITE(iu6,*) "Dimer: Projection"
        WRITE(iu6,*) "Dimer: N*F0", SUM(N*F0)
        WRITE(iu6,*) "Dimer: F0sq", SUM(F0*F0)
      END IF
      ! Project force
      IF (CN<0) THEN
        Feff = F0 - N*SUM(F0*N)*2._q
      ELSE
        Feff = N*SUM(F0*N)*(-1._q)
      END IF
      IF (iu6>=0) THEN
        WRITE(iu6,*) "Dimer: Feffsq ", SUM(Feff*Feff)
        WRITE(iu6,*) "Dimer: Feff"
        WRITE(iu6,'(3ES20.10)') (Feff(:,i), i=1,Nions)
        WRITE(iu6,*) ''
      END IF
    END SUBROUTINE ProjectDimer

! Dimer rotation

    SUBROUTINE RotateDimer()
      REAL(q) a1,a2
      ! Conjugate gradient and modified Newton's method for rotation
      IF (FdStep) THEN
        IF (IU6>=0) WRITE(IU6,*) 'Dimer: Trial Rotation'
        ! Finite difference step
        IF (CGInitFlag) THEN
          CGInitFlag = .FALSE.
          FNold = FN
          GN = FN
        END IF
        a1 = ABS(SUM(FN*FNold))
        a2 = SUM(FNold*FNold)
        IF ((a1<=0.5_q*a2) .AND. (a2/=0._q)) THEN
          GamN = SUM(FN*(FN - FNold))/a2
        ELSE
          IF (IU6>=0) WRITE(IU6,*) 'Dimer: CG reset'
          GamN = 0.0_q
        ENDIF
        IF (IU6>=0) WRITE(IU6,'(A,F14.6)') ' Dimer: Gam',GamN
        GN = FN + GN*GamN
        !GH 02/02/18: need to re-orthogonalize to N
        GN = GN - SUM(GN*N)*N
        !GH 02/02/18
        IF(IU6>=0) WRITE(IU6,'(A,F14.6)') ' Dimer: in RotateDimer, mode norm ',Sqrt(Sum(N*N))
        IF(IU6>=0) WRITE(IU6,'(A,F14.6)') ' Dimer: in RotateDimer N*GN',Sum(GN*N)
        GNu = GN/SQRT(SUM(GN*GN))
        FN1 = SUM(FN*GNu)
        IF (IU6>=0) WRITE(IU6,'(A,F14.6)') ' Dimer: in RotateDimer, N*GNu before rotate',Sum(N*GNu)
        CALL Rotate(N, GNu, PI/4._q)
        IF (IU6>=0) WRITE(IU6,'(A,F14.6)') ' Dimer: in RotateDimer, N*GNu after rotate',Sum(N*GNu)
        R = R0 + N*dR
        FdStep = .FALSE.
        IF (IU6>=0) THEN
          WRITE(IU6,'(A,F14.6,/)') ' Dimer: FN1',FN1
        END IF
        ! Save variables for printing
        CN1 = CN
        FN1r = FNr
      ELSE
        IF (IU6>=0) WRITE(IU6,*) 'Dimer: Rotation'
        RotNum = RotNum + 1
        ! Rotation step
        FN2 = SUM(FN*GNu)
        IF (FN2/=0._q) THEN
          Th = ATAN(FN1/FN2)/(-2._q)
        ELSE
          Th = PI/(-2._q)
        END IF
        IF (FN2>0._q) Th = Th + PI/2._q
        CALL Rotate(N, GNu, Th - PI/4._q)
        CALL SetUnit(N)
        R = R0 + N*dR
        FdStep = .TRUE.
        IF (IU6>=0) THEN
          WRITE(IU6,'(A,F14.6)') ' Dimer: FN2', FN2
          WRITE(IU6,'(A,F14.6,/)') ' Dimer: Th ', Th*180._q/PI
          WRITE(IUdim,'(I5,5F14.5)') Itr, F0r, FN1r, U0, CN1, Th*180._q/PI
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
      V = MOD(V + 100.5_q, 1._q) - 0.5_q
      CALL DIRKAR(Nions,V,A)
    END SUBROUTINE SetPBC

! Rotates both V1 and V2

    SUBROUTINE Rotate(V1,V2,tTh)
      REAL(q),DIMENSION(3,Nions) :: V1,V2,V1tmp
      REAL(q) :: tTh,cTh,sTh
      cTh = COS(tTh)
      sTh = SIN(tTh)
      V1tmp = V1
      V1 = V1*cth + V2*sTh
      V2 = V2*cth - V1tmp*sTh
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
        DO i = 1, Nions
          Vtmp(1:3) = V(1:3,i)
          CALL KARDIR(1, Vtmp, B)
          DO j = 1, 3
            IF (.NOT.DINFO%LSFOR(j,i)) Vtmp(j) = 0._q
          END DO
          CALL DIRKAR(1, Vtmp, A)
          V(1:3,i) = Vtmp(1:3)
        END DO
      ELSE
        Vtmp = SUM(V, DIM=2)/Nions
        DO i = 1,Nions
          V(1:3,i) = V(1:3,i) - Vtmp(1:3)
        END DO
        IF (IU6>=0) WRITE(IU6,'(A7,A10,3F14.6)') 'Dimer:', 'Drift', Vtmp(1:3)
      END IF
    END SUBROUTINE SetConstraints

!**********************************************************************
! IO routines
!**********************************************************************

! Write lowest mode

    SUBROUTINE WriteMode()
      IF (IU6>=0) THEN
        OPEN(210, FILE='NEWMODECAR', ACTION='write', STATUS='replace')
        WRITE(210,'(3ES20.10)') (N(:,i) , i=1,Nions)
        CLOSE(210)
      END IF
    END SUBROUTINE WriteMode

! Read Dimer Variables from the INCAR file

    SUBROUTINE ReadDimerVar()
      INTEGER :: IDUM, IERR, INint
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM

      Seed = 0
      CALL RDATAB(.TRUE.,'INCAR',IU5,'DSeed','=','#',';','I', &
     &            Seed,RDUM,CDUM,LDUM,CHARAC,INint,1,IERR)
      RotMax = 4
      CALL RDATAB(.TRUE.,'INCAR',IU5,'DRotMax','=','#',';','I', &
     &            RotMax,RDUM,CDUM,LDUM,CHARAC,INint,1,IERR)
      dR = 5E-3_q
      CALL RDATAB(.TRUE.,'INCAR',IU5,'DdR','=','#',';','F', &
     &            IDUM,dR,CDUM,LDUM,CHARAC,INint,1,IERR)
      FNMax = 1._q
      CALL RDATAB(.TRUE.,'INCAR',IU5,'DFNMax','=','#',';','F', &
     &            IDUM,FNMax,CDUM,LDUM,CHARAC,INint,1,IERR)
      FNMin = 0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IU5,'DFNMin','=','#',';','F', &
     &            IDUM,FNMin,CDUM,LDUM,CHARAC,INint,1,IERR)

      IF (IU6>=0) THEN
        WRITE(IU6,'(/,A)')     ' Dimer: -----------------'
        WRITE(IU6,'(A)')       ' Dimer: Input Parameters'
        WRITE(IU6,'(A,I7)')    ' Dimer:   Seed',Seed
        WRITE(IU6,'(A,I7)')    ' Dimer: RotMax',RotMax
        WRITE(IU6,'(A,F14.6)') ' Dimer:     dR',dR
        WRITE(IU6,'(A,F14.6)') ' Dimer:  FNMax',FNMax
        WRITE(IU6,'(A,F14.6)') ' Dimer:  FNMin',FNMin
        WRITE(IU6,'(A)')       ' Dimer: -----------------'
      END IF

    END SUBROUTINE ReadDimerVar

! Write the CENTCAR.  This routine is taken from poscar.F with a small change:
! the selective dynamics line is writen to the file.

    SUBROUTINE OUTCENT(IU, LLONG, SZNAM, T_INFO, SCALE, Atmp, NTYP, NITYP, LSDYN, NIONS, POSION, LSFOR)
!      USE prec
!      IMPLICIT REAL(q) (A-H,O-Z)

      LOGICAL LLONG,LSDYN
      TYPE (type_info) :: T_INFO
      REAL(q) :: SCALE
      REAL(q),DIMENSION(3,3) :: Atmp
      REAL(q),DIMENSION(3,NIONS) :: POSION
      LOGICAL,DIMENSION(3,NIONS) :: LSFOR
      INTEGER,DIMENSION(NTYP) :: NITYP
      CHARACTER*40 FORM
      CHARACTER*40 SZNAM
      INTEGER IU, NIONS, NTYP
      INTEGER NT, NI

!-----direct lattice
      WRITE(IU,'(A40)') SZNAM

      WRITE(IU,*)  SCALE
      IF (LLONG) THEN
        FORM = '(1X,3F22.16)'
      ELSE
        FORM = '(1X,3F12.6)'
      ENDIF
      WRITE(IU,FORM) (Atmp(1,I)/SCALE, Atmp(2,I)/SCALE, Atmp(3,I)/SCALE, I=1,3)

      IF (T_INFO%TYPE(1)/='  ') THEN
         WRITE(IU,'(20A5)') (T_INFO%TYPE(NT), NT=1,T_INFO%NTYP)
      ENDIF

      WRITE(IU,'(20I4)') (NITYP(NT), NT=1,NTYP)
      IF (LSDYN) WRITE(IU,'(A18)') 'Selective dynamics'
      WRITE(IU,'(A6)') 'Direct'

      IF (LSDYN) THEN
      IF (LLONG) THEN
        FORM = '(3F20.16,3L4)'
      ELSE
        FORM = '(3F10.6,3L2)'
      ENDIF
      ELSE
      IF (LLONG) THEN
        FORM = '(3F20.16)'
      ELSE
        FORM = '(3F10.6)'
      ENDIF
      ENDIF

      IF (LSDYN) THEN
          WRITE(IU,FORM) &
     &      (POSION(1,NI), POSION(2,NI), POSION(3,NI), &
     &       LSFOR(1,NI), LSFOR(2,NI), LSFOR(3,NI), NI=1,NIONS)
         ELSE
          WRITE(IU,FORM) &
     &      (POSION(1,NI), POSION(2,NI), POSION(3,NI), NI=1,NIONS)
      ENDIF
      IF (.NOT.LLONG) WRITE(IU,*)
      RETURN
    END SUBROUTINE

    END MODULE dimer

