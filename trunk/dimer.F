#include "symbol.inc"
#define DIMER
!**********************************************************************
! RCS:  $Id: dimer.F,v 1.14 2006-10-11 18:22:16 graeme Exp $
!
! This MODULE implements the dimer method.  For more information refer
! to the article J.Chem.Phys 111, 7010 (1999) and the web page
! http://theory.cm.utexas.edu/vtsttools/dimer/
!
!  History
! 
! v0.10, Oct 200 the earliest version
!
! v0.90, May 2001 the first official beta
! - Added some comments 
! - Fixed an inconsistancy in the definition of DMaxNR
! - Fixed the WRITE statements so the code will run on many nodes (this was 
!   a major problem on some platforms) 
! - Added the VASP WFORCE calls so that FILEs are writen more often during 
!   the calculation 
!
! v0.91, July 2001
!
! - Thomas Mattsson <tkmatts@sandia.gov> suggested a modification to the 
!   quickmin (damped dynamics) routine in which any element of the velocity
!   is set to zero if it is misalligned with the force. In some quick tests,
!   this improves the confergence of quickmin by a significant amount.
! - Kiril Tsemekhman <kirilt@ikazki01.chem.washington.edu> fixed a bug so 
!   that the dimer can now be run on many processors. The change was setting
!   the IU6>=0 condition in each RDATAB call to .TRUE. so that every processor
!   reads the dimer parameters, not just the designated IO processor for the image. 
!
! v0.92, June 2002
! - Removed a ' charactor in line 330 which caused some compilers to complain. 
!
! v0.93, Aug 2002
! - Du Yaojun <dyj@mps.ohio-state.edu> pointed out that the vector GNu is not
!   saved between vasp force calculations.  The consequence of this bug is that
!   the minimization direction was off by the finite difference rotation in some
!   iterations.  I have not checked how much this bug was affecting confergence.
!
! v1.00, Aug 2002
! - Incorporated into the vasp source as dimer.F.  The three modules: dimer, neb
!   and dynamical matrix are now separate modules controlled by the chain.F FILE.
!
! v1.01, Mar 2005
! - Formatted the output in dimer.dat and added the headings to out.dat
!
! TO DO:
! - add a single image version of the dimer method
! - add a minimum rotational force variable
!
! Graeme Henkelman
! henkelman@mail.utexas.edu
!
!**********************************************************************
!
    MODULE dimer
      USE base
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
      TYPE(type_info) :: DINFO
      INTEGER :: NM,Nions,Node,IU0,IU6,IUtmp,IUout,IUdim
      INTEGER :: I,NI,NJ,NT,MaxNR,NR
      INTEGER :: ImageNum,Seed,Itr,NNFLAG
      REAL(q),ALLOCATABLE :: Fs(:,:,:),Rs(:,:,:),Vs(:,:,:),Us(:,:)
      REAL(q),ALLOCATABLE :: F(:,:),V(:,:),N(:,:),M(:,:),dV(:,:),N0(:,:)
      REAL(q),ALLOCATABLE :: R(:,:),R0(:,:),R1(:,:),R2(:,:),Gold(:,:)
      REAL(q),ALLOCATABLE :: GN(:,:),G(:,:),Gu(:,:),GNu(:,:),GNold(:,:)
      REAL(q),ALLOCATABLE :: FN(:,:),FNP1(:,:),FNP2(:,:),FNold(:,:)
      REAL(q),ALLOCATABLE :: SPF(:,:),FR(:,:),FN1(:,:),FN2(:,:),FM(:,:)
      REAL(q),ALLOCATABLE :: F1(:,:),F2(:,:),FP1(:,:),FP2(:,:),Fold(:,:)
      REAL(q),ALLOCATABLE :: F1tmp(:,:),F2tmp(:,:),Ntmp(:,:),dRtmp(:,:)
      REAL(q),ALLOCATABLE :: R0tmp(:,:),Gtmp(:,:),GNtmp(:,:),GNutmp(:,:)
      REAL(q),ALLOCATABLE :: Mode(:,:,:)
      REAL(q) :: A(3,3),B(3,3)
      REAL(q) :: U1tmp, U2tmp
      REAL(q) :: Uref,U,U0,U1,U2,CN1,CN,Cth,CR
      REAL(q) :: FD1,FD2,FRr,VRr,FACT
      REAL(q) :: NdT,VdV
      REAL(q) :: NdR,dTh,SimR,FMax,dT,dR,RdR
      REAL(q) :: FNr,FNrp,FNrp1,FNrp2,GamN,Th
      REAL(q) :: Frp,Frp1,Frp2,Gam,MaxMove,MaxFN
      LOGICAL :: CGFlag,NCGFlag,Done,DimerFormed
      LOGICAL :: dThStep,ThStep,dRStep,RdRStep
      LOGICAL :: FirstTranslation,FirstRotation,IonStep,TestFlag
      LOGICAL :: DragFlag,NoDrIFtFlag

    CONTAINS
!
!
!**********************************************************************
!
!  Dimer method:  moves a pair of configurations to a saddle point
!
!**********************************************************************
!
    SUBROUTINE Dimer_Force(optflag,TINFO_NIONS,POSION,TOTEN,TIFOR, & 
                                       LATT_A,LATT_B,IO_IU6)
      LOGICAL :: optflag
      INTEGER :: TINFO_NIONS
      INTEGER :: IO_IU6
      REAL(q) :: POSION(3,NIONS)
      REAL(q) :: TIFOR(3,NIONS),TOTEN
      REAL(q) :: LATT_A(3,3),LATT_B(3,3)

      IF(images==0) RETURN
      ! IF(spring==-1000) RETURN
#if defined(MPI) || defined(MPI_CHAIN)
      U=TOTEN
      F=TIFOR
      R=POSION
      A=LATT_A
      B=LATT_B

      Itr=Itr+1
      ! Convert the position into Cartesian coordinates
      CALL DIRKAR(NIONS,R,A)

      IF(IU6>=0) WRITE(IUdim,'(A20,I7)') 'Iteration',Itr

!     Communicate the between the two images 
      CALL GetFRU()

!     If the dimer is at the first step, it has to be properly
!     formed, ie, have the right spacing (NdR) between the images
      IF(FirstRotation) THEN
        CALL DimerInit()
        IF(.NOT.DimerFormed) THEN
!         Now the dimer is formed, so recalculate forces 
          POSION=R
          CALL KARDIR(NIONS,POSION,B)
          RETURN
        ENDIF
        IF(NCGFlag)THEN
!         Using conjugate gradients for dimer rotation  
          dThStep=.TRUE.
          ThStep=.FALSE.
        ELSE
!         Using quickmin for dimer rotation
          dThStep=.FALSE.
          ThStep=.TRUE.
        ENDIF
!       First rotate the dimer ... translate it later
        dRStep=.FALSE.
        RdRStep=.FALSE.
      ENDIF
      CALL UpDateDimer()
      IF(IU6>=0) CALL WriteDimerData(IUdim)
      IF(IU6>=0) CALL WriteDimerData(IU6)
!     Dimer rotation
      IF(dThStep.OR.ThStep) CALL FindN()
!     Dimer translation
      IF(dRStep.OR.RdRStep) CALL DimerStep()
!     SEND the current position and force back to vasp
      POSION=R
      TIFOR=FR
!     Convert back to direct coordinates (only for positions, not forces)
      CALL KARDIR(NIONS,POSION,B)
      IF(IU6>=0) WRITE(IUdim,*)
!     Empty the FILE buffers
      IF(IU6>=0) CALL WFORCE(IUdim)
      IF(IU6>=0) CALL WFORCE(IUout)
!      CALL flush_(IUdim)
!      CALL flush_(IUout)
#endif
      END SUBROUTINE Dimer_Force
!
!
!**********************************************************************
!
! initialize the dimer
!
!**********************************************************************
!
      SUBROUTINE Dimer_Init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(TYPE_info) :: T_INFO
      EXTERNAL RMARIN

      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt) :: LATT_CUR
      TYPE (TYPE_info) :: T_I
      TYPE (dynamics) :: DYN
      INTEGER IERR,idir

      IF(images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)
      DIO=IO
      DINFO=T_INFO
      IU6=DIO%IU6
      IU0=DIO%IU0
      Nions=T_info%nions
      Node=comm_chain%node_me
!      dT=DYN%POTIM
!     This mass is no good and should be taken out
!      Mass=DINFO%POMASS(1) 
      CALL MAKE_DIR_APP(Node)
      IUtmp=43
      IUout=IUtmp+1
      IUdim=IUtmp+2
      IF(IU6>=0) OPEN(UNIT=IUout, FILE=DIR_APP(1:DIR_LEN)//'out.dat')
      IF(IU6>=0) WRITE(IUout,'(6A16)') 'FMax','FNr','U','CN','Cth','Th'
      IF(IU6>=0) OPEN(UNIT=IUdim, FILE=DIR_APP(1:DIR_LEN)//'dimer.dat')
      Itr=0
      ALLOCATE(Rs(3,Nions,1:images),Vs(3,Nions,1:images))
      ALLOCATE(Fs(3,Nions,1:images),Us(1,1:images))
      ALLOCATE(F(3,Nions),F1(3,Nions),F2(3,Nions),FR(3,Nions))
      ALLOCATE(R(3,Nions),R1(3,Nions),R2(3,Nions),R0(3,Nions))
      ALLOCATE(N(3,Nions),M(3,Nions),SPF(3,Nions),V(3,Nions))
      ALLOCATE(FN(3,Nions),FN1(3,Nions),FN2(3,Nions),dV(3,Nions))
      ALLOCATE(G(3,Nions),Gu(3,Nions),GN(3,Nions),GNu(3,Nions))
      ALLOCATE(FM(3,Nions),FNP1(3,Nions),FNP2(3,Nions))
      ALLOCATE(FNold(3,Nions),GNold(3,Nions),N0(3,Nions))
      ALLOCATE(Fold(3,Nions),Gold(3,Nions))
      ALLOCATE(FP1(3,Nions),FP2(3,Nions))
      ALLOCATE(F1tmp(3,Nions),F2tmp(3,Nions),Ntmp(3,Nions))
      ALLOCATE(Gtmp(3,Nions),GNtmp(3,Nions),GNutmp(3,Nions))
      ALLOCATE(R0tmp(3,Nions),dRtmp(3,Nions))
      IF(Seed/=0) CALL RMARIN(Seed,Seed)
      FirstRotation=.TRUE.
      FirstTranslation=.TRUE.
      Done=.FALSE.
      CALL ReadDimerVar()
      IF(images==1) THEN
        IF(.NOT.DragFlag) THEN
          IF(IU6>=0) WRITE(IU0,*) 'For Images=1 USE Drag.'
          stop
        ENDIF
      ELSEIF(images/=2) THEN
        IF(IU6>=0) WRITE(IU0,*) 'Images must be 2 for the dimer.'
        IF(IU6>=0) WRITE(IU0,*) 'Images must be 1 for Drag.'
        stop
      ENDIF
      IF(DragFlag) THEN
        idir=0
        CALL MAKE_DIR_APP(idir)
        CALL RD_POSCAR_HEAD(LATT_CUR,T_I,NIOND,NIONPD,NTYPD,NTYPPD,IU0,IU6)
        CALL RD_POSCAR(LATT_CUR,T_I,DYN,NIOND,NIONPD,NTYPD,NTYPPD,IU0,IU6)
        IF(T_I%NIONS/=T_INFO%NIONS) THEN
          IF(IU6>=0) THEN
            WRITE(IU0,*) 'ERROR: wrong number of ions in POSCAR FILE'
            stop
          ENDIF
        ENDIF
        R0=DYN%POSION
        IF(DragFlag) R0=MOD(R0+100.5_q,1._q)-0.5_q
        CALL DIRKAR(Nions,R0,LATT_CUR%A)
        node=COMM_CHAIN%NODE_ME
        CALL MAKE_DIR_APP(node)
        IF(DragFlag) N0=RetUnit(R0)
      ENDIF
#endif
!MPI
      END SUBROUTINE Dimer_Init
!
!
!**********************************************************************
!
! Dimer Functions
!
!**********************************************************************
!
! Convert the forces and positions to cartesian coords
! and to reasonable units
!
      SUBROUTINE GetFRU()
      CALL CommR()
      CALL CommF()
      CALL CommU()
      END SUBROUTINE GetFRU
!
! Check IF the dimer is properly formed
!
      FUNCTION DimerCheck()
      LOGICAL :: DimerCheck
      REAL(q) :: R12(3,Nions)
      DimerCheck=.FALSE.
      DimerCheck=(ABS(Dist(R1,R2)-NdR*2._q)<(NdR/1000._q))
      IF(IU6>=0) WRITE(IUdim,'(A20,L7)') 'Check formed',DimerCheck
      IF(DimerCheck)THEN
        N=R1-R2
        CALL SetNoDrIFt(N)
        CALL SetUnit(N)
        CALL SetPN()
      ENDIF
      END FUNCTION DimerCheck
!
! Update the forces and energy of the dimer
!
      SUBROUTINE UpDateDimer()
      CALL SetFN()
      U0=(U1+U2)/2._q+(FD1-FD2)*(NdR/4._q)
      CN=(FD2-FD1)/(2.0*NdR)
      END SUBROUTINE UpDateDimer
!
! Project the forces within the dimer
!
      SUBROUTINE SetFN()
      CALL SetNoDrIFt(F1)
      CALL SetNoDrIFt(F2)
      FD1=SUM(F1*N)
      FD2=SUM(F2*N)
      FN1=F1-N*SUM(F1*N)
      FN2=F2-N*SUM(F2*N)
      FN=(FN1-FN2)/NdR
      FNr=SQRT(SUM(FN**2))
      END SUBROUTINE SetFN
!
! Generate R1 and R2 based upon N
!
      SUBROUTINE SetPN()
      IF(node==1) R=(R0+N*NdR)
      IF(node==2) R=(R0-N*NdR)
      CALL CommR()
      END SUBROUTINE SetPN
!
! Dimer rotation
!
      SUBROUTINE FindN()
      REAL(q) :: FNoldN(3,Nions),FNoldG(3,Nions)
      IF(NCGFlag) THEN
!       Conjugate Gradient and modified Newtons method for rotation
        IF(dThStep) THEN
          IF(FirstRotation) THEN
            GamN=0
            FNold=FN
            GNold=FN
! Used for rotating the dimer in a circle (for testing only)
!            GN=FN
!            GNu=RetUnit(GN)
            FirstRotation=.FALSE.
          ELSE
!           GamN=SUM(FN*(FN-FNold))/SUM(FNold**2)
            GamN=0
          ENDIF
          GN=FN+GNold*GamN
!          GN=RetNoDrIFt(GN)
          GNu=RetUnit(GN)
          FNP1=GNu*SUM(FN*GNu)
          FNrP1=SUM(FNP1*GNu)
          CALL SaveTmp()
          CALL Rotate(N,GNu,dTh)
          IF(IU6>=0) THEN
            WRITE(IUdim,'(A10,A10,F14.6)') 'CG FindN:','dTh',dTh
            WRITE(IUdim,'(A20,F14.6)') 'FNrP1',FNrp1
          ENDIF
!          IF(IU6>=0) WRITE(IUout,*) (Itr-1._q)*dTh,(U1+U2)/2._q,CN,FNrP1,FNr
          CALL SetPN()
          dThStep=.FALSE.
          ThStep=.TRUE.
        ELSE
          NR=NR+1
          FNold=FN
          FNp2=GNu*SUM(FN*GNu)
          FNrp2=SUM(FNp2*GNu)
          Cth=(FNrp1-FNrp2)/dTh
          FNrp=(FNrp1+FNrp2)/2._q
          Th=atan((FNrp/Cth)*2._q)/2._q+dTh/2._q
          IF(Cth<0) Th=Th+PI/2._q
          IF(IU6>=0) THEN
            WRITE(IUdim,'(A10,A10,F14.6)')  'CG FindN:','Th',Th
            WRITE(IUdim,'(A20,2F14.6)') 'FNrP1,2',FNrp1,FNrp2
            WRITE(IUdim,'(A20,F14.6)')  'FNrP',FNrp
            WRITE(IUdim,'(A20,F14.6,/)')  'Cth',Cth
            WRITE(IUout,'(6F16.7)') FMax,FNr,U,CN,Cth,Th
          ENDIF
          CALL RestoreTmp()
          FNoldN=N*SUM(FNold*N)
          FNoldG=GNu*SUM(FNold*GNu)
          FNold=FNold-FNoldN
          FNold=FNold-FNoldG
          CALL Rotate(N,GNu,Th)
!          N=RetUnit(RetNoDrIFt(N))
          N=RetUnit(N)
          CALL SetPN()
          GNold=GNu*SQRT(SUM(GN**2))
          FNold=FNold+N*SQRT(SUM(FNoldN**2))
          FNold=FNold+GNu*SQRT(SUM(FNoldG**2))
          IF((FNr<MaxFN).or.(NR>=MaxNR))THEN
            NR=0
            dRStep=.TRUE.
            ThStep=.FALSE.
          ELSE
            ThStep=.FALSE.
            dThStep=.TRUE.
          ENDIF
        ENDIF
      ELSE
!       QuickMin for dimer rotation
        NR=NR+1
        IF(FirstRotation)THEN
          FirstRotation=.FALSE.
          V=0._q
        ENDIF
        dV=FN*(NdR*NdT)
        VdV=SUM(V*dV)
        IF(IU6>=0) THEN
          WRITE(IUdim,'(A10,A10,F14.6)') 'QM FindN:','Vr(before)',SUM(V*V)
          WRITE(IUdim,'(A20,F14.6)') 'dV^2',SUM(dV*dV)
          WRITE(IUdim,'(A20,F14.6)') 'VdV/SUM(dV**2)',VdV/SUM(dV**2)
        ENDIF
        IF(VdV>0._q) THEN
          V=dV*(1._q+VdV/SUM(dV**2))
        ELSE
          V=dV
        ENDIF
        R1=R1+V*NdT
!        N=RetUnit(RetNoDrIFt(R1-R))
!        N=RetUnit(R1-R0)
        N=R1-R0
        CALL SetPBC(N)
        CALL SetUnit(N)
        CALL SetPN()
        IF(IU6>=0) WRITE(IUdim,'(A20,F14.6)') 'VdV',VdV
        IF(IU6>=0) WRITE(IUout,'(6F16.7)') FMax,FNr,U,CN,Cth,Th
        IF((FNr<MaxFN).or.(NR>=MaxNR))THEN
          NR=0
          dRStep=.TRUE.
          ThStep=.FALSE.
        ENDIF
      ENDIF
!      IF((FNr<MaxFN).or.(NR>=MaxNR))THEN
!        NR=0
!        dRStep=.TRUE.
!        IF(NCGFlag)THEN
!          dThStep=.FALSE.
!        ELSE 
!          ThStep=.FALSE.
!        ENDIF
!      ENDIF
      END SUBROUTINE FindN
!
! Translate the dimer
!
      SUBROUTINE DimerStep()
      IF(CN<0)THEN
         SPF=FR-N*SUM(FR*N)*2._q
      ELSE
        SPF=N*SUM(FR*N)*(-1._q)
!        SPF=FR-N*SUM(FR*N)*2._q
      ENDIF
      IF(CGFlag)THEN
!       Conjugate gradient method for dimer translation
        IF(dRStep)THEN
          dRStep=.FALSE.
          RdRStep=.TRUE.
          IF(FirstTranslation)THEN
            FirstTranslation=.FALSE.
            FirstRotation=.FALSE.
            Gam=0
            Fold=SPF
            Gold=SPF
          ELSE
            Gam=SUM(SPF*(SPF-Fold))/SUM(Fold**2)
          ENDIF
          G=SPF+Gold*Gam
          Gu=RetUnit(G)
          Fold=SPF
          Gold=G
          FP1=Gu*SUM(SPF*Gu)
          FrP1=SUM(FP1*Gu)
          IF(IU6>=0) THEN
            WRITE(IUdim,'(A10,A10,F14.6)') 'dR Step:','dR',dR
            WRITE(IUdim,'(A20,2F14.6)') 'F1,2*Gu',SUM(F1*Gu),SUM(F2*Gu)
            WRITE(IUdim,'(A20,F14.6)') '|SPF|',SQRT(SUM(SPF**2))
            WRITE(IUdim,'(A20,F14.6)') 'Max(SPF)',FMax
          ENDIF
          R0=R0+Gu*dR
        ELSE
          RdRStep=.FALSE.
          IF(NCGFlag) THEN
            dThStep=.TRUE.
          ELSE
            ThStep=.TRUE.
          ENDIF
          FP2=Gu*SUM(SPF*Gu)
          FrP2=SUM(FP2*Gu)
          CR=(FrP1-FrP2)/dR
          FrP=(FrP1+FrP2)/2._q
          IF(CR<0)THEN
            RdR=MaxMove
          ELSE
            RdR=FrP/CR
            IF(ABS(RdR)>MaxMove) RdR=SIGN(MaxMove,RdR)
            RdR=RdR-dR/2._q
          ENDIF
          R0=R0+Gu*RdR
          IF(IU6>=0) THEN
            WRITE(IUdim,'(A10,A10,F14.6)') 'RdR Step:','RdR',RdR
            WRITE(IUdim,'(A20,2F14.6)') 'FrP1,2',FrP1,FrP2
            WRITE(IUdim,'(A20,F14.6)') 'FR*Gu',SUM(FR*Gu)
            WRITE(IUdim,'(A20,F14.6)') 'FrP',FrP
            WRITE(IUdim,'(A20,F14.6)') 'CR',CR
            WRITE(IUdim,'(A20,F14.6)') 'FrP/CR',FrP/CR
          ENDIF
        ENDIF
      ELSE
!       QuickMin for dimer translation
        dRStep=.FALSE.
        IF(NCGFlag) THEN
          dThStep=.TRUE.
        ELSE
          ThStep=.TRUE.
        ENDIF
        dV=SPF*dT
!
!        This implementation zeroes any component of the velocity which
!        is in the opposite direction of the force.  The modIFication
!        which improves the efficiency and stability of QuickMin was
!        suggested by Thomas Mattsson <tkmatts@sandia.gov>.
!
        where(dV*V<0) V=0
!
        VdV=SUM(dV*V)
!        IF(VdV>0._q) THEN
          V=dV*(1._q+VdV/SUM(dV**2))
!        ELSE 
!          V=dV
!        ENDIF
        RdR=SQRT(SUM((V*dT)**2))
        IF(RdR>MaxMove) THEN
          R0=R0+RetUnit(V)*MaxMove
        ELSE
          R0=R0+V*dT
        ENDIF
        IF(IU6>=0) THEN
          WRITE(IUdim,'(A10,A10,F14.6)') 'QM Step:','V',SQRT(SUM(V*V))
          WRITE(IUdim,'(A20,F14.6)') 'dV',SQRT(SUM(dV*dV))
          WRITE(IUdim,'(A20,F14.6)') 'VdV/SUM(dV**2)',VdV/SUM(dV**2)
        ENDIF
      ENDIF
      CALL SetPN()
      END SUBROUTINE DimerStep
!
! Save temporary dimer parameters
!
      SUBROUTINE SaveTmp()
      R0tmp=R0
      Ntmp=N
      F1tmp=F1
      F2tmp=F2
      Gtmp=G
      GNtmp=GN
      GNutmp=GNu
      U1tmp=U1
      U2tmp=U2
      END SUBROUTINE SaveTmp
!
! Restore temporary dimer parameters
!
      SUBROUTINE RestoreTmp()
      R0=R0tmp
      N=Ntmp
      F1=F1tmp
      F2=F2tmp
      G=Gtmp
      GN=GNtmp
      GNu=Gnutmp
      U1=U1tmp
      U2=U2tmp
      CALL SetPN()
      CALL UpDateDimer()
      END SUBROUTINE RestoreTmp
!
!**********************************************************************
!
! Drag Method
!
!**********************************************************************
! 
      SUBROUTINE DragStep()
      IF(IU6>=0) WRITE(IUdim,'(A10,A10,I7)') 'DRAG:','Iteration',Itr
      SPF=F-N0*SUM(F*N0)*2._q
      END SUBROUTINE DragStep
!
!
!**********************************************************************
!
! Communication Functions
!
!**********************************************************************
!
      SUBROUTINE CommF()
      CALL SetSelZero(F)
      Fs(:,:,1:images)=0
      Fs(:,:,node)=F
      CALLMPI_C(M_SUM_d(comm_chain,Fs(1,1,1),nions*3*images))
      F1=Fs(:,:,1)
      F2=Fs(:,:,2)
      FR=(F1+F2)/2._q
      FMax=MAXVAL(ABS(FR))
      END SUBROUTINE CommF

      SUBROUTINE CommV()
      Vs(:,:,1:images)=0
      Vs(:,:,node)=V
      CALLMPI_C(M_SUM_d(comm_chain,Vs(1,1,1),nions*3*images))
      END SUBROUTINE CommV

      SUBROUTINE CommR()
      Rs(:,:,1:images)=0
      Rs(:,:,node)=R
!      IF(IU6>=0) WRITE(IU6,*) 'Node: ',node
      CALLMPI_C(M_SUM_d(comm_chain,Rs(1,1,1),nions*3*images))
      R1=Rs(:,:,1)
      R2=Rs(:,:,2)
!      R0=(R1+R2)/2._q
! GH: fixed this for boundary conditions
      dRtmp=R2-R1
      CALL SetPBC(dRtmp)
      R0=R1+dRtmp/2._q
      END SUBROUTINE CommR

      SUBROUTINE CommU()
      Us(:,1:images)=0
      Us(:,node)=U
      CALLMPI_C(M_SUM_d(comm_chain,Us(1,1),images))
      U1=Us(1,1)
      U2=Us(1,2)
      END SUBROUTINE CommU
!     
!
!**********************************************************************
!
! Vector Functions
!
!**********************************************************************
!
! Finds the distance between R1 and R2
!
      FUNCTION Dist(V1,V2)
      REAL(q) :: V1(3,Nions),V2(3,Nions),V12(3,Nions),Dist
      V12=V2-V1
      CALL KARDIR(Nions,V12,B)
      V12=MOD(V12+100.5_q,1._q)-0.5_q
      CALL DIRKAR(Nions,V12,A)
      Dist=SQRT(SUM(V12**2))
      END FUNCTION Dist
!
! Sets a vector to have the smallest length consistent the the periodic boundaries
! This should REALly be changed to USE the Wigner-Sitz cell
!
      SUBROUTINE SetPBC(V1)
      REAL(q) :: V1(3,Nions)
      CALL KARDIR(Nions,V1,B)
      V1=MOD(V1+100.5_q,1._q)-0.5_q
      CALL DIRKAR(Nions,V1,A)
      END SUBROUTINE SetPBC
!
! Rotates both V1 and V2
!
      SUBROUTINE Rotate(V1,V2,Th)
      REAL(q) :: V1(3,Nions),V2(3,Nions),Th
      REAL(q) :: V1tmp(3,Nions),cth,sth
      cth=COS(Th)
      sth=SIN(Th)
      V1tmp=V1
      V1=V1*cth+V2*sth
      V2=V2*cth-V1tmp*sth
      END SUBROUTINE Rotate
!
! Returns a UNIT vector along V1
!
      FUNCTION RetUnit(V1)
      REAL(q) :: V1(3,Nions)
      REAL(q),dimension(3,Nions) :: RetUnit
      RetUnit=V1*(1._q/SQRT(SUM(V1*V1)))
      END FUNCTION RetUnit
!
! Sets V1 to be a UNIT vector
!
      SUBROUTINE SetUnit(V1)
      REAL(q) :: V1(3,Nions)
      V1=RetUnit(V1)
      END SUBROUTINE SetUnit
!
! Sets V1 to be a random vector
!
      SUBROUTINE SetRand(V1)
      REAL(q) :: V1(3,Nions),RANG
      EXTERNAL RANG
      DO I=1,3
        DO NI=1,NIONS
          V1(I,NI)=RANG(0._q,1._q)
        ENDDO
      ENDDO
      END SUBROUTINE SetRand
!
! Set the constained coordinates in a vector to zero
!
      SUBROUTINE SetSelZero(V1)
      REAL(q) :: V1(3,Nions)
      REAL(q) VTMP(3)
      IF(DINFO%LSDYN) THEN
        DO NI=1,NIONS
          DO NJ=1,3
            VTMP(NJ)=V1(NJ,NI)
          ENDDO
          CALL KARDIR(1,VTMP,B)
          DO NJ=1,3
            IF (.NOT.DINFO%LSFOR(NJ,NI)) VTMP(NJ)=0._q
          ENDDO
          CALL DIRKAR(1,VTMP,A)
          DO NJ=1,3
            V1(NJ,NI)=VTMP(NJ)
          ENDDO
        ENDDO
      ENDIF
      END SUBROUTINE SetSelZero
!
! Remove DrIFt
!
      SUBROUTINE SetNoDrift(V1)
      REAL(q) :: V1(3,Nions)
      REAL(q) VTMP(3)
      IF((.NOT.DINFO%LSDYN).AND.NoDriftFlag) THEN
        VTMP=SUM(V1,DIM=2)
        DO NI=1,NIONS
          V1(:,NI)=V1(:,NI)-VTMP
        ENDDO
        IF(IU6>=0) WRITE(IUdim,'(A20,3F14.6)') 'Drift',VTMP(:)
      ENDIF
      END SUBROUTINE SetNoDrift
!
! Remove Drift
!
      FUNCTION RetNoDrift(V1)
      REAL(q)::V1(3,Nions)
      REAL(q),DIMENSION(3,Nions)::RetNoDrift
      REAL(q) VTMP(3)
      RetNoDrift=V1
      IF(.NOT.DINFO%LSDYN) THEN
        VTMP=SUM(V1,DIM=2)
        DO NI=1,NIONS
          RetNoDrift(:,NI)=V1(:,NI)-VTMP
        ENDDO
        IF(IU6>=0) WRITE(IUdim,'(A20,3F14.6)') 'Drift ',VTMP(:)
      ENDIF
      END FUNCTION RetNoDrift
!
!
!**********************************************************************
!
! IO routines
!
!**********************************************************************
!
! Write initial dimer variables
!
      SUBROUTINE WriteDimerVar(IU)
      INTEGER IU
      WRITE(IU,'(/,A22,/)') 'Dimer Input Parameters'
      WRITE(IU,'(A20,L7)') 'CGFlag',CGFlag
      WRITE(IU,'(A20,L7)') 'NCGFlag',NCGFlag
      WRITE(IU,'(A20,L7)') 'DragFlag',DragFlag
      WRITE(IU,'(A20,I7)') 'Seed',Seed
      WRITE(IU,'(A20,I7)') 'NM',NM
      WRITE(IU,'(A20,I7)') 'MaxNR',MaxNR
      WRITE(IU,'(A20,F14.6)') 'NdR',NdR
      WRITE(IU,'(A20,F14.6)') 'dR',dR
      WRITE(IU,'(A20,F14.6)') 'dth',dth
      WRITE(IU,'(A20,F14.6)') 'SimR',SimR
      WRITE(IU,'(A20,F14.6)') 'MaxMove',MaxMove
      WRITE(IU,'(A20,F14.6)') 'NdT',NdT
      WRITE(IU,'(A20,F14.6)') 'dT',dT
      WRITE(IU,*)
      END SUBROUTINE WriteDimerVar
!
! Write dimer parameters
!
      SUBROUTINE WriteDimerData(IU)
      INTEGER IU
      IF(IU.eq.6) THEN
        WRITE(IU,'(A7,A13,F14.6)') 'Dimer:','CN',CN
        WRITE(IU,'(A7,A13,3F14.6)')'Dimer:','U0,1,2',U0,U1,U2
        WRITE(IU,'(A7,A13,F14.6)') 'Dimer:','|F|',SQRT(SUM(FR**2))
        WRITE(IU,'(A7,A13,2F14.6)')'Dimer:','|F1,2|',SQRT(SUM(F1**2)),SQRT(SUM(F2**2))
        WRITE(IU,'(A7,A13,2F14.6)')'Dimer:','FD1,2',FD1,FD2
        WRITE(IU,'(A7,A13,F14.6)') 'Dimer:','FNr',FNr
        WRITE(IU,'(A7,A13,2F14.6,/)') 'Dimer:','|FN1,2|',SQRT(SUM(FN1**2)),SQRT(SUM(FN2**2))
        WRITE(IU,*) 
      ELSE
        WRITE(IU,'(A20,F14.6)') 'CN',CN
        WRITE(IU,'(A20,3F14.6)') 'U0,1,2',U0,U1,U2
        WRITE(IU,'(A20,F14.6)') '|F|',SQRT(SUM(FR**2))
        WRITE(IU,'(A20,2F14.6)') '|F1,2|',SQRT(SUM(F1**2)),SQRT(SUM(F2**2))
        WRITE(IU,'(A20,2F14.6)') 'FD1,2',FD1,FD2
        WRITE(IU,'(A20,F14.6)') 'FNr',FNr
        WRITE(IU,'(A20,2F14.6,/)') '|FN1,2|',SQRT(SUM(FN1**2)),SQRT(SUM(FN2**2))
      ENDIF
      END SUBROUTINE WriteDimerData
!
!
!**********************************************************************
!
! Dimer initialization FUNCTIONs
!
!**********************************************************************
!
      SUBROUTINE DimerInit()
      ImageNum=0
      IF(IU6>=0) WRITE(IUdim,'(A20)') 'Dimer Init'
      IF(IU6>=0) WRITE(IU6,'(A7,A13)') 'Dimer:','Init'
      DimerFormed=DimerCheck()
      IF(.NOT.DimerFormed) THEN
        IF(IU6>=0) WRITE(IUdim,'(A20,/)') 'Dimer not formed'
        IF(IU6>=0) WRITE(IU6,'(A7,A13,/)') 'Dimer:','not formed'
        IF(ABS(Dist(R1,R2))<(NdR/10._q)) THEN
          CALL SetRand(N)
        ELSE
          IF(IU6>=0) THEN
            WRITE(IUdim,'(A20)') 'Forming from images'
            WRITE(IU6,'(A7,A20)') 'Dimer:','forming from images'
          ENDIF
          N=R1-R2
        ENDIF
!        CALL SetNoDrift(N)
        CALL SetPBC(N)
        CALL SetSelZero(N)
        CALL SetUnit(N)
        CALL CommR()
        CALL SetPN()
      ENDIF
      END SUBROUTINE DimerInit
!
! Read all Dimer Variables from the INCAR FILE
! -- and write them to the OUTCAR FILE
!
      SUBROUTINE ReadDimerVar()
      INTEGER :: IDUM,IERR,Nint
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM
      IF(IU6>=0) WRITE(IUdim,'(A20)') 'Loading variables'

      CGFlag=.TRUE.
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DCGFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,CGFlag,CHARAC,Nint,1,IERR)
      NCGFlag=.TRUE.
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DNCGFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,NCGFlag,CHARAC,Nint,1,IERR)
      NoDrIFtFlag=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'NoDrIFtFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,NoDrIFtFlag,CHARAC,Nint,1,IERR)
      DragFlag=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DragFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,DragFlag,CHARAC,Nint,1,IERR)
      NM=1
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DNM','=','#',';','I', &
     &            NM,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      Seed=0
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DSeed','=','#',';','I', &
     &            Seed,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      MaxNR=1
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DMaxNR','=','#',';','I', &
     &            MaxNR,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      NdR=1E-3_q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DNdR','=','#',';','F', &
     &            IDUM,NdR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dR=1E-3_q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DdR','=','#',';','F', &
     &            IDUM,dR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dTh=2E-3_q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DdTh','=','#',';','F', &
     &            IDUM,dTh,CDUM,LDUM,CHARAC,Nint,1,IERR)
      SimR=1E-2_q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DSimR','=','#',';','F', &
     &            IDUM,SimR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      MaxMove=0.1_q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DMaxMove','=','#',';','F', &
     &            IDUM,MaxMove,CDUM,LDUM,CHARAC,Nint,1,IERR)
      NdT=1._q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DNdT','=','#',';','F', &
     &            IDUM,NdT,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dT=0.1_q
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DdT','=','#',';','F', &
     &            IDUM,dT,CDUM,LDUM,CHARAC,Nint,1,IERR)
      MaxFN=1
      CALL RDATAB(.TRUE.,'INCAR',DIO%IU5,'DMaxFN','=','#',';','I', &
     &            IDUM,MaxFN,CDUM,LDUM,CHARAC,Nint,1,IERR)
      IonStep=.FALSE.
      NR=0
      IF(IU6>=0) CALL WriteDimerVar(IUdim)
      IF(IU6>=0) CALL WriteDimerVar(IU6)
      END SUBROUTINE ReadDimerVar
!
! Load modes from a FILE (not currently USEd)
!
      SUBROUTINE LoadModes()
      REAL(q) :: V1(3,Nions)
      CHARACTER(15) :: ImageFileName,ImageNumText
      ImageNum=ImageNum+1
      IF(IU6>=0) WRITE(ImageNumText,*) ImageNum
      ImageFileName = "Image"//TRIM(ImageNumText(2:))//".dat"
      OPEN(UNIT=IUtmp, FILE='images/'//ImageFileName)
      IF(IU6>=0) WRITE(IUtmp,*) (V1(:,ni),ni=1,nions)
      IF(IU6>=0) CLOSE(UNIT=IUtmp)
      END SUBROUTINE LoadModes

      END MODULE dimer
