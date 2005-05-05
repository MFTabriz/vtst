#include "symbol.inc"
#define DIMER
!**********************************************************************
! RCS:  $Id: dimer.F,v 1.1.1.1 2005-05-05 22:30:51 graeme Exp $
!
! This module implements the dimer method.  For more information refer
! to the article J.Chem.Phys 111, 7010 (1999) and the web page
! http://ikazki01.chem.washington.edu/vasp/index.html
!
! This is the first beta release of this method.  There will be some bugs.
! The code was developed on an IBM AIX system.  Some problems were found 
! porting this code to other platforms.  Let me know if you have problems,
! or questions.
!
!  History
! 
! v0.10, Oct 200 the earliest version
!
!
! v0.90, May 2001 the first official beta
! - Added some comments 
! - Fixed an inconsistancy in the definition of DMaxNR
! - Fixed the write statements so the code will run on many nodes (this was 
!   a major problem on some platforms) 
! - Added the VASP WFORCE calls so that files are writen more often during 
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
!   the IU6>=0 condition in each RDATAB call to .True. so that every processor
!   reads the dimer parameters, not just the designated IO processor for the image. 
!
! v0.92, June 2002
! - Removed a quote charactor in line 330 which caused some compilers to complain. 
!
! v0.93, Aug 2002
! - Du Yaojun <dyj@mps.ohio-state.edu> pointed out that the vector GNu is not
!   saved between vasp force calculations.  The consequence of this bug is that
!   the minimization direction was off by the finite difference rotation in some
!   iterations.  I have not checked how much this bug was affecting confergence.
!
! v1.00, Aug 2002
! - Incorporated into the vasp source as dimer.F.  The three modules: dimer, neb
!   and dynamical matrix are now separate modules controlled by the chain.F file.
!
! v1.01, May 2004
! - Added a steepest descent reset rule to the conjugate gradient routine
!   as suggested by Andri Arnaldsson
!
! Graeme Henkelman
! graeme@u.washington.edu
!
!**********************************************************************
!
      module dimer
      use base
      use prec                                                              
      use main_mpi
      use poscar
      use lattice
      use constant
      implicit none
      save
      private
      public :: Dimer_Force,Dimer_Init

      type(in_struct) :: DIO
      type(type_info) :: DINFO
      integer :: NM,Nions,Node,IU0,IU6,IUtmp,IUout,IUdim
      integer :: I,NI,NJ,NT,MaxNR,NR
      integer :: ImageNum,Seed,Itr,TestItr,NNFLAG
      real(q),allocatable :: Fs(:,:,:),Rs(:,:,:),Vs(:,:,:),Us(:,:)
      real(q),allocatable :: F(:,:),V(:,:),N(:,:),M(:,:),dV(:,:),N0(:,:)
      real(q),allocatable :: R(:,:),R0(:,:),R1(:,:),R2(:,:),Gold(:,:)
      real(q),allocatable :: GN(:,:),G(:,:),Gu(:,:),GNu(:,:),GNold(:,:)
      real(q),allocatable :: FN(:,:),FNP1(:,:),FNP2(:,:),FNold(:,:)
      real(q),allocatable :: SPF(:,:),FR(:,:),FN1(:,:),FN2(:,:),FM(:,:)
      real(q),allocatable :: F1(:,:),F2(:,:),FP1(:,:),FP2(:,:),Fold(:,:)
      real(q),allocatable :: F1tmp(:,:),F2tmp(:,:),Ntmp(:,:)
      real(q),allocatable :: R0tmp(:,:),Gtmp(:,:),GNtmp(:,:),GNutmp(:,:)
      real(q),allocatable :: Mode(:,:,:)
!      real(q),allocatable :: Ream(:),Feam(:)
      real(q) :: A(3,3),B(3,3)
!      real(q) :: spring=10
      real(q) :: U1tmp, U2tmp, a1, a2
      real(q) :: Uref,U,U0,U1,U2,CN1,CN,Cth,CR
      real(q) :: FD1,FD2,FRr,VRr,FACT
      real(q) :: NdT,VdV,Alpha
      real(q) :: NdR,dTh,SimR,FMax,dT,dR,RdR
      real(q) :: FNr,FNrp,FNrp1,FNrp2,GamN,Th
      real(q) :: Frp,Frp1,Frp2,Gam,MaxMove,MaxFN
      logical :: ContFlag,CGFlag,NCGFlag,Done,DimerFormed
      logical :: dThStep,ThStep,dRStep,RdRStep
      logical :: FirstTranslation,FirstRotation,IonStep,TestFlag
      logical :: ImageFlag,ARTFlag,DragFlag,NoDriftFlag,EAMFlag
!
!
!**********************************************************************
!
!  routine for forces between the images on the elastic band
!
!**********************************************************************
!
      contains
!      subroutine DIMER_FORCE(T_INFO,DYN,LATT_CUR,TOTEN,TIFOR)
!      subroutine CHAIN_FORCE(T_INFO%NIONS,DYN%POSION,TOTEN,TIFOR,LATT_CUR%A,LATT_CUR%B,IO%IU6)
!      subroutine CHAIN_FORCE(TINFO_NIONS,POSION,TOTEN,TIFOR,LATT_A,LATT_B,IO_IU6)
      subroutine DIMER_FORCE(TINFO_NIONS,POSION,TOTEN,TIFOR,LATT_A,LATT_B,IO_IU6)
      integer :: TINFO_NIONS
      integer :: IO_IU6
      real(q) :: POSION(3,NIONS)
      real(q) :: TIFOR(3,NIONS),TOTEN
      real(q) :: LATT_A(3,3),LATT_B(3,3)

!      external FORCE
      if(images==0) return
!      if(spring==-1000) return
#if defined(MPI) || defined(MPI_CHAIN)
      U=TOTEN
      F=TIFOR
      R=POSION
      A=LATT_A
      B=LATT_B

      Itr=Itr+1
!     Convert the position into Cartesian coordinates
      call DIRKAR(NIONS,R,A)

!     This is only used for the ART or DRAG method
!     and not with the dimer method
      if(ARTFlag.or.DragFlag)then
        if(ARTFlag) call ARTStep()
        if(DRAGFlag) call DragStep()
        TIFOR=SPF
        if(ImageFlag) then
          if(IU6>=0) call WriteImage(R)
!         if(IU6>=0) call WriteNormal(R0tmp)
!         if(IU6>=0) call WriteForce(FR)
        endif
        return
      endif	   

      if(IU6>=0) write(IUdim,*) 'Dimer: Iteration ',Itr
!      if(IU6>=0) write(IU0,*) 'Dimer: A ',A(1,1),A(2,2),A(3,3)
!      if(EAMFlag) call EAMForce()

!     Communicate the between the two images 
      call GetFRU()

!     If the dimer is at the first step, it has to be properly
!     formed, ie, have the right spacing (NdR) between the images
      if(FirstRotation) then
        call DimerInit()
        if(.not.DimerFormed) then
!         Now the dimer is formed, so recalculate forces 
          POSION=R
          call KARDIR(NIONS,POSION,B)
          return
        endif
        if(NCGFlag)then
!         Using conjugate gradients for dimer rotation  
          dThStep=.True.
          ThStep=.False.
        else
!         Using quickmin for dimer rotation
          dThStep=.False.
          ThStep=.True.
        endif
!       First rotate the dimer ... translate it later
        dRStep=.False.
        RdRStep=.False.
      endif
      call UpDateDimer()
      if(IU6>=0) call WriteDimerData(IUdim)
      if(IU6>=0) call WriteDimerData(IU6)
!     Dimer rotation
      if(dThStep.or.ThStep) call FindN()
!     Dimer translation
      if(dRStep.or.RdRStep) call DimerStep()
!      This was an idea to let vasp minimize the position of the dimer
!      It ended up being too complicated so it is not used
!      if(IonStep)then
!        TIFOR=SPF
!        R=R0
!        if(Node==1.and.ImageFlag) then
!          Only generate images if i am the first node
!          call WriteImage(R0)
!        endif
!      endif         
!     Send the current position and force back to vasp
      POSION=R
	TIFOR=FR
!     Convert back to direct coordinates (only for positions, not forces)
      call KARDIR(NIONS,POSION,B)
      if(IU6>=0) write(IUdim,*)
!     Empty the file buffers
      if(IU6>=0) call WFORCE(IUdim)
      if(IU6>=0) call WFORCE(IUout)
!      call flush_(IUdim)
!      call flush_(IUout)
#endif
      end subroutine Dimer_Force
!
!
!**********************************************************************
!
! initialize the chain (repeated image mode)
! read the spring constant
! and  the two outer images, these images are kept fixed
! during the entire simulation
!
!**********************************************************************
!
!      subroutine Chain_Init(T_INFO,IO)
      subroutine Dimer_Init(T_INFO,IO)
      type(in_struct) :: IO
      type(type_info) :: T_INFO
      external RMARIN

      integer NIOND,NIONPD,NTYPPD,NTYPD
      type (latt) :: LATT_CUR
      type (type_info) :: T_I
      type (dynamics) :: DYN
      integer IERR,idir

      if(images==0) return
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
      call MAKE_DIR_APP(Node)
      IUtmp=43
      IUout=IUtmp+1
      IUdim=IUtmp+2
      if(IU6>=0) open(unit=IUout, file=DIR_APP(1:DIR_LEN)//'out.dat')
      if(IU6>=0) open(unit=IUdim, file=DIR_APP(1:DIR_LEN)//'dimer.dat')
      Itr=0
 !     TestItr=0
      allocate(Rs(3,Nions,1:images),Vs(3,Nions,1:images))
      allocate(Fs(3,Nions,1:images),Us(1,1:images))
      allocate(F(3,Nions),F1(3,Nions),F2(3,Nions),FR(3,Nions))
      allocate(R(3,Nions),R1(3,Nions),R2(3,Nions),R0(3,Nions))
      allocate(N(3,Nions),M(3,Nions),SPF(3,Nions),V(3,Nions))
      allocate(FN(3,Nions),FN1(3,Nions),FN2(3,Nions),dV(3,Nions))
      allocate(G(3,Nions),Gu(3,Nions),GN(3,Nions),GNu(3,Nions))
      allocate(FM(3,Nions),FNP1(3,Nions),FNP2(3,Nions))
      allocate(FNold(3,Nions),GNold(3,Nions),N0(3,Nions))
      allocate(Fold(3,Nions),Gold(3,Nions))
      allocate(FP1(3,Nions),FP2(3,Nions))
      allocate(F1tmp(3,Nions),F2tmp(3,Nions),Ntmp(3,Nions))
      allocate(Gtmp(3,Nions),GNtmp(3,Nions),GNutmp(3,Nions))
      allocate(R0tmp(3,Nions))
!      allocate(Feam(3*Nions),Ream(3*Nions))
      if(Seed/=0) call RMARIN(Seed,Seed)
      FirstRotation=.True.
      FirstTranslation=.True.
      Done=.False.
      call ReadDimerVar()
      if(images==1)then
        if(ARTFlag.eqv.DragFlag) then
          if(IU6>=0) write(IU0,*)'For Images=1 set only ARTFlag or DragFlag'
          stop
        endif
      elseif(images/=2)then
        if(IU6>=0) write(IU0,*)'It is required that Images=2 for the dimer method.'
        if(IU6>=0) write(IU0,*)'Or Images=1 for the ART or Drag methods.'
        stop
      endif
      if(ARTFlag.OR.DragFlag)then
        idir=0
        call MAKE_DIR_APP(idir)
        call RD_POSCAR_HEAD(LATT_CUR,T_I,NIOND,NIONPD,NTYPD,NTYPPD,IU0,IU6)
        call RD_POSCAR(LATT_CUR,T_I,DYN,NIOND,NIONPD,NTYPD,NTYPPD,IU0,IU6)
        if(T_I%NIONS/=T_INFO%NIONS) then
          if(IO%IU0>=0) then
            if(IU6>=0) write(IU0,*)'ERROR: image mode number of ions wrong'
            stop
          endif
        endif
        R0=DYN%POSION
        if(DragFlag) R0=mod(R0+100.5_q,1._q)-0.5_q
        call DIRKAR(Nions,R0,LATT_CUR%A)
        node=COMM_CHAIN%NODE_ME
        call MAKE_DIR_APP(node)
        if(DragFlag) N0=RetUnit(R0)
      endif
!      if(EAMFlag) call PotInit(7.1,7.25);
#endif
!MPI
      end subroutine Dimer_Init
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
      subroutine GetFRU()
      call CommR()
      call CommF()
      call CommU()
      end subroutine GetFRU
!
! Check if the dimer is properly formed
!
      function DimerCheck()
      logical :: DimerCheck
      real(q) :: R12(3,Nions)
      DimerCheck=.false.
      DimerCheck=(abs(Dist(R1,R2)-NdR*2._q)<(NdR/1000._q))
      if(DimerCheck)then
        if(IU6>=0) write(IUdim,*) 'Dimer: Dimer Check = True'
        N=R1-R2
        call SetNoDrift(N)
        call SetUnit(N)
        call SetPN()
      else 
        if(IU6>=0) write(IUdim,*) 'Dimer: Dimer Check = False'
      endif
      end function DimerCheck
!
! Update the forces and energy of the dimer
!
      subroutine UpDateDimer()
      call SetFN()
      U0=(U1+U2)/2._q+(FD1-FD2)*(NdR/4._q)
      CN=(FD2-FD1)/(2.0*NdR)
      end subroutine UpDateDimer
!
! Project the forces within the dimer
!
      subroutine SetFN()
      call SetNoDrift(F1)
      call SetNoDrift(F2)
      FD1=sum(F1*N)
      FD2=sum(F2*N)
      FN1=F1-N*sum(F1*N)
      FN2=F2-N*sum(F2*N)
      FN=(FN1-FN2)/NdR
      FNr=sqrt(sum(FN**2))
      end subroutine SetFN
!
! Generate R1 and R2 based upon N
!
      subroutine SetPN()
      if(node==1) R=(R0+N*NdR)
      if(node==2) R=(R0-N*NdR)
      call CommR()
      end subroutine SetPN

!
! Dimer rotation
!
      subroutine FindN()
      real(q) :: FNoldN(3,Nions),FNoldG(3,Nions)
      if(NCGFlag)then
!       Conjugate Gradient and modified Newtons method for rotation
        if(dThStep)then
          if(FirstRotation)then
            GamN=0
            FNold=FN
            GNold=FN
! Used for rotating the dimer in a circle (for testing only)
!            GN=FN
!            GNu=RetUnit(GN)
            FirstRotation=.False.
          else
            a1=abs(sum(FN*FNold)) 
            a2=sum(FNold*FNold)
            if(a1<=0.5*a2)then
              GamN=sum(FN*(FN-FNold))/sum(FNold**2)
            else
              GamN=0
            endif
          endif
          GN=FN+GNold*GamN
!          GN=RetNoDrift(GN)
          GNu=RetUnit(GN)
          FNP1=GNu*sum(FN*GNu)
          FNrP1=sum(FNP1*GNu)
          call SaveTmp()
          call Rotate(N,GNu,dTh)
          if(IU6>=0) write(IUdim,*) 'Dimer: FindN, dTh ',dTh
          if(IU6>=0) write(IUdim,*) 'Dimer: FNrP1 ',FNrp1
!          if(IU6>=0) write(IUout,*) (Itr-1._q)*dTh,(U1+U2)/2._q,CN,FNrP1,FNr
          call SetPN()
          dThStep=.False.
          ThStep=.True.
        else
          NR=NR+1
          FNold=FN
          FNp2=GNu*sum(FN*GNu)
          FNrp2=sum(FNp2*GNu)
          Cth=(FNrp1-FNrp2)/dTh
          FNrp=(FNrp1+FNrp2)/2._q
          Th=atan((FNrp/Cth)*2._q)/2._q+dTh/2._q
          if(Cth<0) Th=Th+PI/2._q
          if(IU6>=0) write(IUdim,*) 'Dimer: FindN, Th ',Th
          if(IU6>=0) write(IUdim,*) 'Dimer: FNrP1,2 ',FNrp1,FNrp2
          if(IU6>=0) write(IUdim,*) 'Dimer: FNrP ',FNrp
          if(IU6>=0) write(IUdim,*) 'Dimer: Cth ',Cth
          if(IU6>=0) write(IUout,'(6E15.7)') FMax,FNr,U,CN,Cth,Th
          call RestoreTmp()
          FNoldN=N*sum(FNold*N)
          FNoldG=GNu*sum(FNold*GNu)
          FNold=FNold-FNoldN
          FNold=FNold-FNoldG
          call Rotate(N,GNu,Th)
!          N=RetUnit(RetNoDrift(N))
          N=RetUnit(N)
          call SetPN()
          GNold=GNu*sqrt(sum(GN**2))
          FNold=FNold+N*sqrt(sum(FNoldN**2))
          FNold=FNold+GNu*sqrt(sum(FNoldG**2))
          if((FNr<MaxFN).or.(NR>=MaxNR))then
            NR=0
            dRStep=.True.
            ThStep=.False.
          else
            ThStep=.False.
	      dThStep=.True.
          endif
        endif
      else
!       QuickMin for dimer rotation
        NR=NR+1
        if(FirstRotation)then
          FirstRotation=.False.
          V=0._q
        endif
        dV=FN*(NdR*NdT)
        VdV=sum(V*dV)
        if(IU6>=0) write(IUdim,*) 'Dimer: FindN, QM '
        if(IU6>=0) write(IUdim,*) 'Dimer: Vr(before) ',sum(V*V)
        if(IU6>=0) write(IUdim,*) 'Dimer: dV^2 ',sum(dV*dV)
        if(IU6>=0) write(IUdim,*) 'Dimer: VdV/sum(dV**2) ',VdV/sum(dV**2)
        if(VdV>0._q) then
          V=dV*(1._q+VdV/sum(dV**2))
        else
          V=dV
        endif
        R1=R1+V*NdT
!        N=RetUnit(RetNoDrift(R1-R))
        N=RetUnit(R1-R0)
        call SetPN()
        if(IU6>=0) write(IUdim,*) 'Dimer: VdV ',VdV
        if(IU6>=0) write(IUout,'(6E15.7)') FMax,FNr,U,CN,Cth,Th
        if((FNr<MaxFN).or.(NR>=MaxNR))then
          NR=0
          dRStep=.True.
          ThStep=.False.
        endif
      endif
      end subroutine FindN
!
! Translate the dimer
!
      subroutine DimerStep()
      if(CN<0)then
         SPF=FR-N*sum(FR*N)*2._q
      else
        SPF=N*sum(FR*N)*(-1._q)
!        SPF=FR-N*sum(FR*N)*2._q
      endif
      if(CGFlag)then
!       Conjugate gradient method for dimer translation
        if(dRStep)then
          dRStep=.False.
          RdRStep=.True.
          if(FirstTranslation)then
            FirstTranslation=.False.
            FirstRotation=.False.
            Gam=0
            Fold=SPF
            Gold=SPF
          else
            a1=abs(sum(SPF*Fold))
            a2=sum(Fold*Fold)
            if(a1<=0.5*a2)then
              Gam=sum(SPF*(SPF-Fold))/sum(Fold**2)
            else
              GamN=0
            endif
          endif
          G=SPF+Gold*Gam
          Gu=RetUnit(G)
          Fold=SPF
          Gold=G
          FP1=Gu*sum(SPF*Gu)
          FrP1=sum(FP1*Gu)
          if(IU6>=0) write(IUdim,*) 'Dimer: Step dR   ',dR
          if(IU6>=0) write(IUdim,*) 'Dimer: F1,2*Gu   ',sum(F1*Gu),sum(F2*Gu)
          if(IU6>=0) write(IUdim,*) 'Dimer: |SPF|     ',sqrt(sum(SPF**2))
          if(IU6>=0) write(IUdim,*) 'Dimer: Max(SPF)  ',FMax
          R0=R0+Gu*dR
        else
          RdRStep=.False.
          if(NCGFlag) then
            dThStep=.True.
          else
            ThStep=.True.
          endif
          FP2=Gu*sum(SPF*Gu)
          FrP2=sum(FP2*Gu)
          CR=(FrP1-FrP2)/dR
          FrP=(FrP1+FrP2)/2._q
          if(CR<0)then
            RdR=MaxMove
          else
            RdR=FrP/CR
            if(abs(RdR)>MaxMove) RdR=Sign(MaxMove,RdR)
            RdR=RdR-dR/2._q
          endif
          R0=R0+Gu*RdR
          if(IU6>=0) write(IUdim,*) 'Dimer: Step, RdR ',RdR
          if(IU6>=0) write(IUdim,*) 'Dimer: FrP1,2    ',FrP1,FrP2
          if(IU6>=0) write(IUdim,*) 'Dimer: FR*Gu     ',sum(FR*Gu)
          if(IU6>=0) write(IUdim,*) 'Dimer: FrP       ',FrP
          if(IU6>=0) write(IUdim,*) 'Dimer: CR        ',CR
          if(IU6>=0) write(IUdim,*) 'Dimer: FrP/CR    ',FrP/CR
          if(Node==1.and.ImageFlag) then
            if(IU6>=0) call WriteImage(R0)
!           if(IU6>=0) call WriteNormal(N)
!           if(IU6>=0) call WriteForce(FR)
          endif
        endif
      else
!       QuickMin for dimer translation
        dRStep=.False.
        if(NCGFlag) then
          dThStep=.True.
        else
          ThStep=.True.
        endif
        dV=SPF*dT
!
!        This implementation zeroes any component of the velocity which
!        is in the opposite direction of the force.  The modification
!        which improves the efficiency and stability of QuickMin was
!        suggested by Thomas Mattsson <tkmatts@sandia.gov>.
!
        where(dV*V<0) V=0
!
        VdV=sum(dV*V)
!        if(VdV>0._q) then
          V=dV*(1._q+VdV/sum(dV**2))
!        else 
!          V=dV
!        endif
        RdR=sqrt(sum((V*dT)**2))
        if(RdR>MaxMove) then
          R0=R0+RetUnit(V)*MaxMove
        else
          R0=R0+V*dT
        endif
        if(IU6>=0) write(IUdim,*) 'Dimer: Step, QM '
        if(IU6>=0) write(IUdim,*) 'Dimer: V ',sqrt(sum(V*V))
        if(IU6>=0) write(IUdim,*) 'Dimer: dV ',sqrt(sum(dV*dV))
        if(IU6>=0) write(IUdim,*) 'Dimer: VdV/sum(dV**2) ',VdV/sum(dV**2)
        if(Node==1.and.ImageFlag) then
          if(IU6>=0) call WriteImage(R0)
!          if(IU6>=0) call WriteNormal(N)
!          if(IU6>=0) call WriteForce(FR)
        endif
      endif
      call SetPN()
      end subroutine DimerStep
!
! Save temporary dimer parameters
!
      subroutine SaveTmp()
      R0tmp=R0; Ntmp=N
      F1tmp=F1; F2tmp=F2
      Gtmp=G; GNtmp=GN; GNutmp=GNu;
      U1tmp=U1; U2tmp=U2
      end subroutine SaveTmp
!
! Restore temporary dimer parameters
!
      subroutine RestoreTmp()
      R0=R0tmp; N=Ntmp
      F1=F1tmp; F2=F2tmp
      G=Gtmp; GN=GNtmp; GNu=Gnutmp
      U1=U1tmp; U2=U2tmp
      call SetPN()
      call UpDateDimer()
      end subroutine RestoreTmp
!
! EAM force evaluation
!
!      subroutine EAMForce()
!      external Force
!      do NI=1,NIONS
!         do NJ=1,3
!            Ream(3*(NI-1)+NJ)=R(NJ,NI)
!         enddo
!      enddo
!      call FORCE(Nions,Ream,Feam,U,A(1,1),A(2,2),A(3,3))
!      do NI=1,NIONS
!         do NJ=1,3
!            F(NJ,NI)=Feam(3*(NI-1)+NJ)
!         enddo
!      enddo
!      end subroutine EAMForce
!
!
!**********************************************************************
!
! Art Method
!
!**********************************************************************
! 
      subroutine ARTStep()
      if(IU6>=0) write(IUdim,*) 'ART: Iteration ',Itr
      R0tmp=R-R0
      call KARDIR(Nions,R0tmp,B)
      R0tmp=mod(R0tmp+100.5_q,1._q)-0.5_q
      call DIRKAR(Nions,R0tmp,A)
      call SetUnit(R0tmp)
      SPF=F-R0tmp*sum(F*R0tmp)*(1._q+Alpha)
      end subroutine ARTStep
!
!
!**********************************************************************
!
! Drag Method
!
!**********************************************************************
! 
      subroutine DragStep()
      if(IU6>=0) write(IUdim,*) 'DRAG: Iteration ',Itr
      SPF=F-N0*sum(F*N0)*2._q
      end subroutine DragStep
!
!
!**********************************************************************
!
! Communication Functions
!
!**********************************************************************
!
      subroutine CommF()
      call SetSelZero(F)
      Fs(:,:,1:images)=0
      Fs(:,:,node)=F
      CALLMPI_C(M_sum_d(comm_chain,Fs(1,1,1),nions*3*images))
      F1=Fs(:,:,1);F2=Fs(:,:,2)
      FR=(F1+F2)/2._q
      FMax=maxval(abs(FR))
      end subroutine CommF

      subroutine CommV()
      Vs(:,:,1:images)=0
      Vs(:,:,node)=V
      CALLMPI_C(M_sum_d(comm_chain,Vs(1,1,1),nions*3*images))
      end subroutine CommV

      subroutine CommR()
      Rs(:,:,1:images)=0
      Rs(:,:,node)=R
      if(IU6>=0) write(IU6,*) 'Node: ',node
      CALLMPI_C(M_sum_d(comm_chain,Rs(1,1,1),nions*3*images))
      R1=Rs(:,:,1);R2=Rs(:,:,2);
      R0=(R1+R2)/2._q
      end subroutine CommR

      subroutine CommU()
      Us(:,1:images)=0
      Us(:,node)=U
      CALLMPI_C(M_sum_d(comm_chain,Us(1,1),images))
      U1=Us(1,1);U2=Us(1,2)
      end subroutine CommU
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
      function Dist(V1,V2)
      real(q) :: V1(3,Nions),V2(3,Nions),V12(3,Nions),Dist
      V12=V2-V1
      call KARDIR(Nions,V12,B)
      V12=mod(V12+100.5_q,1._q)-0.5_q
      call DIRKAR(Nions,V12,A)
      Dist=sqrt(sum(V12**2))
      end function Dist
!
! Rotates both V1 and V2
!
      subroutine Rotate(V1,V2,Th)
      real(q) :: V1(3,Nions),V2(3,Nions),Th
      real(q) :: V1tmp(3,Nions),cth,sth
      cth=cos(Th)
      sth=sin(Th)
      V1tmp=V1
      V1=V1*cth+V2*sth
      V2=V2*cth-V1tmp*sth
      end subroutine Rotate
!
! Returns a unit vector along V1
!
      function RetUnit(V1)
      real(q) :: V1(3,Nions)
      real(q),dimension(3,Nions) :: RetUnit
      RetUnit=V1*(1._q/sqrt(sum(V1*V1)))
      end function RetUnit
!
! Sets V1 to be a unit vector
!
      subroutine SetUnit(V1)
      real(q) :: V1(3,Nions)
      V1=RetUnit(V1)
      end subroutine SetUnit
!
! Sets V1 to be a random vector
!
      subroutine SetRand(V1)
      real(q) :: V1(3,Nions),RANE
	external RANE
      do I=1,3
        do NI=1,NIONS
          V1(I,NI)=RANE()
        enddo
      enddo
      end subroutine SetRand
!
! Set the constained coordinates in a vector to zero
!
      subroutine SetSelZero(V1)
      real(q) :: V1(3,Nions)
      real(q) VTMP(3)
      if(DINFO%LSDYN) then
        do NI=1,NIONS
          do NJ=1,3
            VTMP(NJ)=V1(NJ,NI)
          enddo
          call KARDIR(1,VTMP,B)
          do NJ=1,3
            if (.NOT.DINFO%LSFOR(NJ,NI)) VTMP(NJ)=0._q
          enddo
          call DIRKAR(1,VTMP,A)
          do NJ=1,3
            V1(NJ,NI)=VTMP(NJ)
          enddo
        enddo
      endif
      end subroutine SetSelZero
!
! Remove Drift
!
      subroutine SetNoDrift(V1)
      real(q) :: V1(3,Nions)
      real(q) VTMP(3)
      if((.NOT.DINFO%LSDYN).AND.NoDriftFlag) then
        VTMP=sum(V1,DIM=2)
        do NI=1,NIONS
          V1(:,NI)=V1(:,NI)-VTMP
        enddo
        if(IU6>=0) write(IUdim,*) 'Dimer: Drift ',VTMP(:)
      endif
      end subroutine SetNoDrift
!
! Remove Drift
!
      function RetNoDrift(V1)
      real(q)::V1(3,Nions)
      real(q),dimension(3,Nions)::RetNoDrift
      real(q) VTMP(3)
      RetNoDrift=V1
      if(.NOT.DINFO%LSDYN) then
        VTMP=sum(V1,DIM=2)
        do NI=1,NIONS
          RetNoDrift(:,NI)=V1(:,NI)-VTMP
        enddo
        if(IU6>=0) write(IUdim,*) 'Dimer: Drift ',VTMP(:)
      endif
      end function RetNoDrift
!
!
!**********************************************************************
!
! IO routines
!
!**********************************************************************
!
! Writes vector out to OUTCAR
!
      subroutine WriteOut(V1,Name,IU)
      real(q) :: V1(3,Nions)
      character(*) :: Name
      integer IU
      write(IU,*) 'Dimer: ',Name
      write(IU,'(3F13.5)') (V1(:,ni),ni=1,nions)
      end subroutine WriteOut
!
! Writes vector out to OUTCAR (in exponential notation)
!
      subroutine WriteOutE(V1,Name,IU)
      real(q) :: V1(3,Nions)
      character(*) :: Name
      integer IU
      write(IU,*) 'Dimer: ',Name
      write(IU,'(3E13.5)') (V1(:,ni),ni=1,nions)
      end subroutine WriteOutE
!
! Writes position to image file (only if node==1)
!
      subroutine WriteImage(V1)
      real(q) :: V1(3,Nions)
      character(15) :: ImageFileName,ImageNumText
      ImageNum=ImageNum+1
      write(ImageNumText,*) ImageNum
      ImageFileName = "Image"//Trim(ImageNumText(2:))//".dat"
      open(unit=IUtmp, file='images/'//ImageFileName)
      write(IUtmp,'(3F13.5)') (V1(:,ni),ni=1,nions)
!      write(IUtmp,*) (V1(:,ni),ni=1,nions)
      close(unit=IUtmp)
      end subroutine WriteImage
!
! Writes normal to normal file (only if node==1)
!
      subroutine WriteNormal(V1)
      real(q) :: V1(3,Nions)
      character(15) :: NormalFileName,ImageNumText
      write(ImageNumText,*) ImageNum
      NormalFileName = "Normal"//Trim(ImageNumText(2:))//".dat"
      open(unit=IUtmp, file='images/'//NormalFileName)
      write(IUtmp,'(3F13.5)') (V1(:,ni),ni=1,nions)
      close(unit=IUtmp)
      end subroutine WriteNormal
!
! Writes force to force file (only if node==1)
!
      subroutine WriteForce(V1)
      real(q) :: V1(3,Nions)
      character(15) :: ForceFileName,ImageNumText
      write(ImageNumText,*) ImageNum
      ForceFileName = "Force"//Trim(ImageNumText(2:))//".dat"
      open(unit=IUtmp, file='images/'//ForceFileName)
      write(IUtmp,'(3F13.5)') (V1(:,ni),ni=1,nions)
      close(unit=IUtmp)
      end subroutine WriteForce
!
! Write initial dimer variables
!
      subroutine WriteDimerVar(IU)
      integer IU
      write(IU,*) 
      write(IU,*) 'Dimer: Input Parameters'
      write(IU,*) 
      write(IU,'(A10,L5)') 'ContFlag',ContFlag
      write(IU,'(A10,L5)') 'CGFlag',CGFlag
      write(IU,'(A10,L5)') 'NCGFlag',NCGFlag
      write(IU,'(A10,L5)') 'IFlag',ImageFlag
      write(IU,'(A10,L5)') 'ARTFlag',ARTFlag
      write(IU,'(A10,L5)') 'DragFlag',DragFlag
      write(IU,'(A10,I5)') 'Seed',Seed
      write(IU,'(A10,I5)') 'NM',NM
      write(IU,'(A10,I5)') 'MaxNR',MaxNR
      write(IU,'(A10,F11.5)') 'NdR',NdR
      write(IU,'(A10,F11.5)') 'dR',dR
      write(IU,'(A10,F11.5)') 'dth',dth
      write(IU,'(A10,F11.5)') 'SimR',SimR
      write(IU,'(A10,F11.5)') 'MaxMove',MaxMove
      write(IU,'(A10,F11.5)') 'NdT',NdT
      write(IU,'(A10,F11.5)') 'dT',dT
!      write(IU,'(A10,F11.5)') 'Mass',Mass
      write(IU,*)
      end subroutine WriteDimerVar
!
! Write dimer parameters
!
      subroutine WriteDimerData(IU)
      integer IU
      write(IU,*) 'Dimer: CN      ',CN
      write(IU,*) 'Dimer: U0,1,2: ',U0,U1,U2
      write(IU,*) 'Dimer: |F|     ',sqrt(sum(FR**2))
      write(IU,*) 'Dimer: |F1,2|  ',sqrt(sum(F1**2)),sqrt(sum(F2**2))
      write(IU,*) 'Dimer: FD1,2   ',FD1,FD2
      write(IU,*) 'Dimer: FNr     ',FNr
      write(IU,*) 'Dimer: |FN1,2| ',sqrt(sum(FN1**2)),sqrt(sum(FN2**2))
      write(IU,*)
      end subroutine WriteDimerData
!
!
!**********************************************************************
!
! Dimer initialization functions
!
!**********************************************************************
!
      subroutine DimerInit()
      ImageNum=0
      if(IU6>=0) write(IUdim,*) 'Dimer: Dimer Init'
      if(IU6>=0) write(IU6,*) 'Dimer: Dimer Init'
      DimerFormed=DimerCheck()
      if(.not.DimerFormed) then
        if(IU6>=0) write(IUdim,*) 'Dimer: Dimer Not Formed'
        if(IU6>=0) write(IUdim,*)
        if(IU6>=0) write(IU6,*) 'Dimer: Dimer Not Formed'
        if(IU6>=0) write(IU6,*)
        if(abs(Dist(R1,R2))<(NdR/10._q)) then
          call SetRand(N)
          R=Rs(:,:,1)
        else
          if(IU6>=0) write(IUdim,*) 'Dimer: Forming dimer from images'
          if(IU6>=0) write(IU6,*) 'Dimer: Forming dimer from images'
          N=R1-R2
        endif
        call SetNoDrift(N)
        call SetSelZero(N)
        call SetUnit(N)
        call CommR()
        call SetPN()
      endif
      end subroutine DimerInit
!
! Read all Dimer Variables from the INCAR file
! -- and write them to the OUTCAR file
!
      subroutine ReadDimerVar()
      integer :: IDUM,IERR,Nint
      character*1 :: CHARAC
      complex(q) :: CDUM 
      logical :: LDUM
      real(q) :: RDUM
      if(IU6>=0) write(IUdim,*) 'Dimer: Loading Dimer Variables'

      ContFlag=.False.
      call RDATAB(.True.,'INCAR',DIO%IU5,'DContFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ContFlag,CHARAC,Nint,1,IERR)
      CGFlag=.True.
      call RDATAB(.True.,'INCAR',DIO%IU5,'DCGFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,CGFlag,CHARAC,Nint,1,IERR)
      NCGFlag=.True.
      call RDATAB(.True.,'INCAR',DIO%IU5,'DNCGFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,NCGFlag,CHARAC,Nint,1,IERR)
      ImageFlag=.False.
      call RDATAB(.True.,'INCAR',DIO%IU5,'DIFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ImageFlag,CHARAC,Nint,1,IERR)
      NoDriftFlag=.False.
      call RDATAB(.True.,'INCAR',DIO%IU5,'NoDriftFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,NoDriftFlag,CHARAC,Nint,1,IERR)
      ARTFlag=.False.
      call RDATAB(.True.,'INCAR',DIO%IU5,'ARTFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ARTFlag,CHARAC,Nint,1,IERR)
      DragFlag=.False.
      call RDATAB(.True.,'INCAR',DIO%IU5,'DragFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,DragFlag,CHARAC,Nint,1,IERR)
      EAMFlag=.False.
      call RDATAB(.True.,'INCAR',DIO%IU5,'DEAMFlag','=','#',';','L', &
     &            IDUM,RDUM,CDUM,EAMFlag,CHARAC,Nint,1,IERR)
      NM=1
      call RDATAB(.True.,'INCAR',DIO%IU5,'DNM','=','#',';','I', &
     &            NM,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      Seed=0
      call RDATAB(.True.,'INCAR',DIO%IU5,'DSeed','=','#',';','I', &
     &            Seed,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      MaxNR=1
      call RDATAB(.True.,'INCAR',DIO%IU5,'DMaxNR','=','#',';','I', &
     &            MaxNR,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      NdR=1E-3_q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DNdR','=','#',';','F', &
     &            IDUM,NdR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dR=1E-3_q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DdR','=','#',';','F', &
     &            IDUM,dR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dTh=2E-3_q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DdTh','=','#',';','F', &
     &            IDUM,dTh,CDUM,LDUM,CHARAC,Nint,1,IERR)
      SimR=1E-2_q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DSimR','=','#',';','F', &
     &            IDUM,SimR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      MaxMove=0.1_q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DMaxMove','=','#',';','F', &
     &            IDUM,MaxMove,CDUM,LDUM,CHARAC,Nint,1,IERR)
      NdT=1._q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DNdT','=','#',';','F', &
     &            IDUM,NdT,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dT=0.1_q
      call RDATAB(.True.,'INCAR',DIO%IU5,'DdT','=','#',';','F', &
     &            IDUM,dT,CDUM,LDUM,CHARAC,Nint,1,IERR)
      MaxFN=1
      call RDATAB(.True.,'INCAR',DIO%IU5,'DMaxFN','=','#',';','I', &
     &            IDUM,MaxFN,CDUM,LDUM,CHARAC,Nint,1,IERR)
      Alpha=1.0
      call RDATAB(.True.,'INCAR',DIO%IU5,'DAlpha','=','#',';','I', &
     &            IDUM,Alpha,CDUM,LDUM,CHARAC,Nint,1,IERR)
      IonStep=.False.
      NR=0
      if(IU6>=0) call WriteDimerVar(IUdim)
      if(IU6>=0) call WriteDimerVar(IU6)
      end subroutine ReadDimerVar
!
! Load modes from a file (not currently used)
!
      subroutine LoadModes()
      real(q) :: V1(3,Nions)
      character(15) :: ImageFileName,ImageNumText
      ImageNum=ImageNum+1
      if(IU6>=0) write(ImageNumText,*) ImageNum
      ImageFileName = "Image"//Trim(ImageNumText(2:))//".dat"
      open(unit=IUtmp, file='images/'//ImageFileName)
      if(IU6>=0) write(IUtmp,*) (V1(:,ni),ni=1,nions)
      if(IU6>=0) close(unit=IUtmp)
      end subroutine LoadModes

      end module dimer
