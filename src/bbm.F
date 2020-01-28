#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: bbm.F,v 1.33 2007/06/06 16:03:09 Liang Zhang Exp $
!
! This module implements the bond boost method.  For more information refer
! to the web page:
!   http://theory.cm.utexas.edu/vtsttools/
! and the articles 
!   [R. A. Miron and K. A. Fichthorn, J. Chem. Phys. 119, 6210 (2003)]
!
! Liang Zhang
! brightzhang@mail.utexas.edu
!
!**********************************************************************

    MODULE bbm
      USE prec
      USE main_mpi
      USE poscar
      USE lattice
      USE constant 

      IMPLICIT NONE
      SAVE 
      PRIVATE
      PUBLIC :: bbm_step,bbm_init

      TYPE(latt) :: L
      INTEGER :: Nions,IU0,IU5,IU6,NBAS,NRMS,NBBS,RMDS,NTBS
      INTEGER :: i,j,n
      INTEGER,ALLOCATABLE,DIMENSION(:) :: BALIST,RMLIST
      REAL(q),ALLOCATABLE,DIMENSION(:) :: AVE_T_BOND,ORI_BOND
      INTEGER,ALLOCATABLE,DIMENSION(:,:) :: T_RATM,ATOMR
      REAL(q) :: QRR,PRR,DVMAX,RCUT,BPOTIM,POTIM,BTEBEG
      LOGICAL :: LBBM

    CONTAINS

!**********************************************************************
! Initialize the bbm
!**********************************************************************

    SUBROUTINE bbm_init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO
    
      INTEGER ::  FLAG,ACOUNT,NI,NJ

      IU0 = IO%IU0
      IU5 = IO%IU5
      IU6 = IO%IU6
      Nions = T_info%nions

      WRITE(IU6,'(/,A)') "====================Bond-Boost Method Initializing===================="
      CALL ReadBBMVar()

      NTBS=NBAS*(NBAS-1)/2+NBAS*(NIONS-NBAS)
      ALLOCATE(AVE_T_BOND(NTBS))
      ALLOCATE(T_RATM(NTBS,2)) 


      AVE_T_BOND=0

      NRMS=NIONS-NBAS
      ALLOCATE(RMLIST(NRMS))
      ACOUNT=0
      DO I=1,Nions
        FLAG=1
        DO J=1,NBAS
          IF(I==BALIST(J)) THEN
            FLAG=0
          ENDIF
        ENDDO
        IF(FLAG == 1) THEN
          ACOUNT=ACOUNT+1
          RMLIST(ACOUNT)=I
        ENDIF
      ENDDO
      
      
      IF (NRMS /= ACOUNT) THEN
         WRITE(IU6,'(A)') "NRMS does not equal to counted number !"
      ENDIF
     
    END SUBROUTINE bbm_Init

!**********************************************************************
!  Bond-Boost Method For each ionic step
!**********************************************************************

    SUBROUTINE bbm_Step(OptFlag,POSION,TOTEN,TIFOR,LATT_A,LATT_B)
      LOGICAL :: OptFlag
      INTEGER :: I,J,k,BCOUNT,RMSTEPS,NSTEP=0,NREG=0,FLAG2=0
      INTEGER :: ATOMR_1,ATOMR_2,BI,BATM1,BATM2,MI,MATOM1,MATOM2
      REAL(q),DIMENSION(7000) :: RBOND
      REAL(q),DIMENSION(NTBS) ::DBOND
      INTEGER,DIMENSION(7000,2) :: ATOMR_TMP
      Real(q),Allocatable,DIMENSION(:) :: EPSR_Q,BOOST_BOND
      REAL(q),DIMENSION(3,NIONS) :: POSION,TIFOR,TADF,ORI_FOR
      REAL(q),DIMENSION(3,3) :: LATT_A,LATT_B
      REAL(q) :: TOTEN,EPSR_MAX,A_EPS_M,BOOST_FACT,AVE_BOOST_FACT,SUM_V
      REAL(q) :: DFORCER,MFORCE,VTMP(3)
      REAL(q) :: SDTIME=0.0,SPTIME=0.0,SDTIME_B=0.0 
      SAVE SDTIME,SPTIME,SDTIME_B,NSTEP,NREG
  
      
            
!     WRITE(IU6,'(/,A,I6,A)') '======================== Running Bond-Boost( NSTEP=',NSTEP,' ) ============================='
       
      SDTIME=SDTIME+BPOTIM
      
      IF(NREG .LE. RMDS) THEN
        FLAG2=0
      ELSE
        FLAG2=1
      ENDIF

      IF (FLAG2 == 0) THEN
        CALL RMDSTEP(SPTIME,NREG,SDTIME,POSION,LATT_A,DBOND)
        AVE_T_BOND=AVE_T_BOND+DBOND/RMDS
      ELSEIF (FLAG2 == 1) THEN
        IF (NREG == RMDS+1) THEN
          BCOUNT=0
          CALL BSTSELECT(RBOND,ATOMR_TMP,BCOUNT)
          NBBS=BCOUNT
           
          ALLOCATE(ORI_BOND(NBBS),ATOMR(NBBS,2))

          DO I=1,NBBS
            ORI_BOND(I)=RBOND(I)
            !WRITE(IU6,*)"NREG=",NREG,"ORI_BOND=",ORI_BOND(I),"RBOND=",RBOND(I)
            ATOMR(I,1)=ATOMR_TMP(I,1)
            ATOMR(I,2)=ATOMR_TMP(I,2)
          ENDDO           
        ENDIF

        ALLOCATE(BOOST_BOND(NBBS),EPSR_Q(NBBS))

        EPSR_Q=0
        BOOST_BOND=0
        ORI_FOR=TIFOR
        NREG=NREG+1
        NSTEP=NSTEP+1
        CALL BOOSTSTEP(NBBS,POSION,TIFOR,LATT_A,ATOMR,BOOST_BOND,EPSR_Q,EPSR_MAX,A_EPS_M,SUM_V,BOOST_FACT,&
          &AVE_BOOST_FACT,SDTIME,SPTIME,SDTIME_B,TADF,ORI_FOR,MI,MATOM1,MATOM2,MFORCE)
        
        IF (LBBM) THEN
          CALL BoostOut(NSTEP,NREG,POSION,TIFOR,BOOST_BOND,EPSR_Q,EPSR_MAX,A_EPS_M,SUM_V,BOOST_FACT,&
            &AVE_BOOST_FACT,SDTIME,SPTIME,SDTIME_B,TADF,ORI_FOR,MI,MATOM1,MATOM2,MFORCE)
        ENDIF

      ENDIF
      WRITE(IU6,'(/,A,I6,A,/)') '===================== End of Bond-Boost Session( NSTEP=',NSTEP,' ) ========================='
    END SUBROUTINE bbm_Step

!==============================================================
! Routine for Boost steps 
!==============================================================
   
    SUBROUTINE BoostStep(NBBS,POSION,TIFOR,LATT_A,ATOMR,BOOST_BOND,EPSR_Q,EPSR_MAX,A_EPS_M,SUM_V,BOOST_FACT,&
                     &AVE_BOOST_FACT,SDTIME,SPTIME,SDTIME_B,TADF,ORI_FOR,MI,MATOM1,MATOM2,MFORCE)

      INTEGER :: I,J,k,BCOUNT,NBBS
      INTEGER :: ATOMR_1,ATOMR_2,MI,MATOM1,MATOM2
      INTEGER,DIMENSION(NBBS,2) :: ATOMR
      Real(q),DIMENSION(NBBS) :: EPSR_Q,BOOST_BOND
      Real(q),DIMENSION(3,NBBS) :: ADD_FOR
      REAL(q),DIMENSION(3,NIONS) :: POSION,TIFOR,TADF,ORI_FOR,FIN_FOR
      REAL(q),DIMENSION(3,3) :: LATT_A
      REAL(q) :: BOND,EPSR_MAX,A_EPS_M,BOOST_FACT,AVE_BOOST_FACT,SUM_V
      REAL(q) :: DFORCER,MFORCE,FACT_1,FACT_2,Q_X,Q_Y,Q_Z,VTMP(3)
      REAL(q) :: SDTIME,SPTIME,SDTIME_B
       
          DO I=1,NBBS
             ATOMR_1= ATOMR(I,1)
             ATOMR_2= ATOMR(I,2)
             CALL BANearest(ATOMR_2,ATOMR_1,POSION,LATT_A,Q_X,Q_Y,Q_Z,BOND)
             BOOST_BOND(I)=BOND
          ENDDO
          
    !EPSR_Q(I): Fractional change in bond length from equilibrium of Boosted Bond(I)
          DO I=1,NBBS
              EPSR_Q(I)=(BOOST_BOND(I)-ORI_BOND(I))/(QRR*ORI_BOND(I))
          ENDDO

    !Find the Max{|Eps(i)|}
          EPSR_MAX=0.0
          DO I=1,NBBS
              EPSR_MAX=MAX(EPSR_MAX,ABS(EPSR_Q(I)))
          ENDDO

    !A_EPS_M: Envelope Function A[EPSR_MAX]         
         A_EPS_M=(1.0_q-EPSR_MAX**2)**2/(1.0_q-PRR**2*EPSR_MAX**2)

         IF(EPSR_MAX .LT. 1.0_q) THEN
             SUM_V=0
             DO I=1,NBBS
              SUM_V=SUM_V+DVMAX*(1.0_q-EPSR_Q(I)**2)/REAL(NBBS)
             ENDDO
             BOOST_FACT=A_EPS_M*SUM_V
             SPTIME=SPTIME+BPOTIM*EXP(BOOST_FACT*1E5_q/BTEBEG/8.63125_q)
             SDTIME_B=SDTIME_B+BPOTIM
             AVE_BOOST_FACT=(SPTIME-(SDTIME-SDTIME_B))/SDTIME_B
              
     
             DO I=1,NIONS
                DO J=1,3
                   TADF(J,I)=0.0_q   !TADF : Total Added Force
                ENDDO
             ENDDO

!Calcualte the Added force by boost potential
             DO I=1,NBBS
                IF(ABS(EPSR_Q(I)) .LT. EPSR_MAX) THEN
                   FACT_1=2.0_q*A_EPS_M*DVMAX*EPSR_Q(I)/QRR/ORI_BOND(I)/REAL(NBBS)
                   DFORCER=FACT_1 !DFORCER=delta_F(I)
                ELSEIF (ABS(EPSR_Q(I)) .EQ. EPSR_MAX) THEN
                   FACT_1=2.0_q*A_EPS_M*DVMAX*EPSR_Q(I)/QRR/ORI_BOND(I)/REAL(NBBS)
                   FACT_2=2.0_q*(1.0_q-(EPSR_Q(I))**2)*EPSR_Q(I)*(2.0_q*(1.0_q-PRR**2*(EPSR_Q(I))**2 )-PRR**2*(1.0_q-(EPSR_Q(I))**2))/QRR/ORI_BOND(I)/(1.0_q-PRR**2*(EPSR_Q(I))**2)**2
                   DFORCER=FACT_1+SUM_V*FACT_2  !DFORCER=delta_F(I)
                ENDIF

                ATOMR_1=ATOMR(I,1)
                ATOMR_2=ATOMR(I,2)
                CALL BANearest(ATOMR_2,ATOMR_1,POSION,LATT_A,Q_X,Q_Y,Q_Z,BOND)
                  
                !Project the Force to x,y,z directions
                ADD_FOR(1,I)=Q_X/BOOST_BOND(I)*DFORCER
                ADD_FOR(2,I)=Q_Y/BOOST_BOND(I)*DFORCER
                ADD_FOR(3,I)=Q_Z/BOOST_BOND(I)*DFORCER

                TADF(1,ATOMR_1)=TADF(1,ATOMR_1)+ADD_FOR(1,I)
                TADF(2,ATOMR_1)=TADF(2,ATOMR_1)+ADD_FOR(2,I)
                TADF(3,ATOMR_1)=TADF(3,ATOMR_1)+ADD_FOR(3,I)

                TADF(1,ATOMR_2)=TADF(1,ATOMR_2)-ADD_FOR(1,I)
                TADF(2,ATOMR_2)=TADF(2,ATOMR_2)-ADD_FOR(2,I)
                TADF(3,ATOMR_2)=TADF(3,ATOMR_2)-ADD_FOR(3,I)
            ENDDO

            ORI_FOR=TIFOR
            FIN_FOR=ORI_FOR+TADF          
            TIFOR=FIN_FOR  

      ELSEIF (EPSR_MAX .GE. 1.0_q) THEN
           BOOST_FACT=0
           SPTIME=SPTIME+BPOTIM
           AVE_BOOST_FACT=(SPTIME-(SDTIME-SDTIME_B))/SDTIME_B
      ENDIF
      WRITE(IU6,'(A,F20.10,3X,A,F20.10,3X,A,F20.10)')"Simulation Time=",SDTIME,"Physical Time=",SPTIME,"Simulation Time in Bond-Boost Steps=",SDTIME_B

   END SUBROUTINE Booststep

!==============================================================
! Read BBM Variables from the INCAR file
!==============================================================

    SUBROUTINE ReadBBMVar()
      INTEGER :: IDUM,IERR,INint,NI
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM 
      LOGICAL :: LDUM
      REAL(q) :: RDUM

     NBAS =0
      CALL RDATAB(.true.,'INCAR',IU5,'NBAS','=','#',';','I', &
     &            NBAS,RDUM,CDUM,LDUM,CHARAC,INint,1,IERR)
     
     ALLOCATE(BALIST(NBAS))
     DO NI=1,NBAS
         BALIST(NI)=0
      ENDDO
     CALL RDATAB(.true.,'INCAR',IU5,'BALIST','=','#',';','I', &
     &            BALIST,RDUM,CDUM,LDUM,CHARAC,INint,NBAS,IERR)
     RMDS =2000
      CALL RDATAB(.true.,'INCAR',IU5,'RMDS','=','#',';','I', &
     &            RMDS,RDUM,CDUM,LDUM,CHARAC,INint,1,IERR)
     BTEBEG=300.0_q
      CALL RDATAB(.true.,'INCAR',IU5,'TEBEG','=','#',';','F', &
     &            IDUM,BTEBEG,CDUM,LDUM,CHARAC,N,1,IERR)
     POTIM=1.0_q
      CALL RDATAB(.true.,'INCAR',IU5,'POTIM','=','#',';','F', &
     &            IDUM,POTIM,CDUM,LDUM,CHARAC,N,1,IERR)
     BPOTIM=POTIM
      CALL RDATAB(.true.,'INCAR',IU5,'BPOTIM','=','#',';','F', &
     &            IDUM,BPOTIM,CDUM,LDUM,CHARAC,N,1,IERR)
     QRR =0.5_q
      CALL RDATAB(.true.,'INCAR',IU5,'QRR','=','#',';','F', &
     &            IDUM,QRR,CDUM,LDUM,CHARAC,INint,1,IERR)
     PRR =0.90_q 
      CALL RDATAB(.true.,'INCAR',IU5,'PRR','=','#',';','F', &
     &            IDUM,PRR,CDUM,LDUM,CHARAC,INint,1,IERR)
     RCUT =3.0_q
      CALL RDATAB(.true.,'INCAR',IU5,'RCUT','=','#',';','F', &
     &            IDUM,RCUT,CDUM,LDUM,CHARAC,INint,1,IERR)
     LBBM =.FALSE.
      CALL RDATAB(.true.,'INCAR',IU5,'LBBM','=','#',';','L', &
     &            IDUM,RDUM,CDUM,LBBM,CHARAC,INint,1,IERR)
     DVMAX =0.3_q 
      CALL RDATAB(.true.,'INCAR',IU5,'DVMAX','=','#',';','F', &
     &            IDUM,DVMAX,CDUM,LDUM,CHARAC,INint,1,IERR)

      IF (IU6>=0) THEN
        WRITE(IU6,'(/,A)')     '  --------------Bond-Boost Input Parameters--------------'
        WRITE(IU6,'(A,11X,I7)')    ' Bond-Boost:  NBAS',NBAS
!        WRITE(IU6,'(A,9X,<NBAS>I5)')   ' Bond-Boost:  BALIST',BALIST(:)
        WRITE(IU6,'(A,9X,999I5)')   ' Bond-Boost:  BALIST',BALIST(:)
        WRITE(IU6,'(A,10X,F14.6)') ' Bond-Boost:  POTIM',BPOTIM
        WRITE(IU6,'(A,10X,F14.6)') ' Bond-Boost:  TEBEG',BTEBEG
        WRITE(IU6,'(A,12X,F14.6)') ' Bond-Boost:  PRR',PRR
        WRITE(IU6,'(A,12X,F14.6)') ' Bond-Boost:  QRR',QRR
        WRITE(IU6,'(A,11X,F14.6)') ' Bond-Boost:  RCUT',RCUT
        WRITE(IU6,'(A,9X,F14.6)') ' Bond-Boost:  DV_MAX',DVMAX 
        WRITE(IU6,'(A,I7)') ' Bond-Boost:  Reguler MD Steps',RMDS
        IF (LBBM) THEN
          WRITE(IU6,'(A)') 'Bond-Boost: Detailed Info will be written'
        ENDIF
      END IF
    
   END SUBROUTINE ReadBBMVar

!==============================================================
! Routine for Regular MD steps to get Equilibrium Bonds lengths
!==============================================================
   
  SUBROUTINE RMDSTEP(SPTIME,NREG,SDTIME,POSION,LATT_A,DBOND)
      INTEGER :: NREG,RMSTEPS
      INTEGER :: ATOMR_1,ATOMR_2
      REAL(q),DIMENSION(NTBS):: DBOND,BOND_TMP
      INTEGER,DIMENSION(NTBS,2) :: ATOMR_TMP
      REAL(q),DIMENSION(3,NIONS) :: POSION
      REAL(q),DIMENSION(3,3) :: LATT_A
      REAL(q) :: SDTIME,SPTIME
      
      SPTIME=SPTIME+BPOTIM
      RMSTEPS=RMDS-NREG
      NREG=NREG+1
  
!      WRITE(IU6,'(A,F10.2,A,I4,A)')'Regular MD has been performed for ', SDTIME,' fs, Bond Boost will start after ',RMSTEPS,' steps'            
!      WRITE(IU6,'(A,/)')'-------------------------------------------------------------------'
      CALL Bond_ini(POSION,LATT_A,BOND_TMP,ATOMR_TMP)
      DBOND=BOND_TMP
      T_RATM=ATOMR_TMP
                 
  END SUBROUTINE RMDSTEP
!=================================================================
!
!Nearest Neighbor between 2 atoms with PBC 
!
!=================================================================
    SUBROUTINE BANearest(SN1,SN2,POSION,L,Q_X,Q_Y,Q_Z,BOND)
    USE prec
    USE lattice
    USE main_mpi
    IMPLICIT NONE
 
    INTEGER IU6,IU0,NBAS,SN1,SN2
    REAL(q),DIMENSION(3,3) :: L
    REAL(q),DIMENSION(3,NIONS) :: POSION
    REAL(q) :: D_x,D_y,D_z,Q_x,Q_y,Q_z,BOND

!Direct coordinate  
    D_x=POSION(1,SN1)-POSION(1,SN2)
    D_y=POSION(2,SN1)-POSION(2,SN2)
    D_z=POSION(3,SN1)-POSION(3,SN2)

!PBC
        IF(D_x .GT. 0.5_q) THEN
            D_x=D_x-1.0_q
        ELSEIF(D_x .LT. -0.5_q) THEN
            D_x=D_x+1.0_q
        ELSE 
            D_x=D_x
        ENDIF

        IF(D_y .GT. 0.5_q) THEN
            D_y=D_y-1.0_q
        ELSEIF(D_y .LT. -0.5_q) THEN
            D_y=D_y+1.0_q
        ELSE
            D_y=D_y
        ENDIF
               
        IF(D_z .GT. 0.5_q) THEN
            D_z=D_z-1.0_q
        ELSEIF(D_z .LT. -0.5_q) THEN
            D_z=D_z+1.0_q
        ELSE 
            D_z=D_z
        ENDIF

!Cartesian Coordinate   
     Q_x=L(1,1)*D_x+L(1,2)*D_y+L(1,3)*D_z
     Q_y=L(2,1)*D_x+L(2,2)*D_y+L(2,3)*D_z
     Q_z=L(3,1)*D_x+L(3,2)*D_y+L(3,3)*D_z
    
     BOND=SQRT(Q_x**2+Q_y**2+Q_z**2)
     !WRITE(IU6,*)"TEST BONDnearst=",BOND
   
    END SUBROUTINE BANearest

!===================================================================
!Record the Bond length and the label of corresponding bonded atoms
!===================================================================     
    SUBROUTINE Bond_ini(POSION,LATT_A,DBOND,ATOMR_TMP)
   
      INTEGER :: I,J,CCOUNT
      REAL(q),DIMENSION(3,NIONS) :: POSION
      REAL(q),DIMENSION(3,3) :: LATT_A
      REAL(q),DIMENSION(NTBS) ::  DBOND
      INTEGER,DIMENSION(NTBS,2) :: ATOMR_TMP
      REAL(q) :: BOND,Q_X,Q_Y,Q_Z
 
!Record bonds between Boosted Atoms   
    CCOUNT=0
!    WRITE(IU6,'(2X,A,7X,A,7X,A,7X,A)')'Bond','Atom1','Atom2','Bond Length'
    DO I=1,NBAS
         DO J=I+1,NBAS
             CALL BANearest(BALIST(I),BALIST(J),POSION,LATT_A,Q_X,Q_Y,Q_Z,BOND)
             CCOUNT=CCOUNT+1
!             WRITE(IU6,'(I4,2I12,F18.6)') CCOUNT,BALIST(I),BALIST(J),BOND
              DBOND(CCOUNT)=BOND
              ATOMR_TMP(CCOUNT,1)=BALIST(I)
              ATOMR_TMP(CCOUNT,2)=BALIST(J)
            ENDDO
      ENDDO

!Record bonds between Boosted Atoms and Rest Atoms
    DO I=1,NBAS
       DO J=1,NRMS
           CALL BANearest(BALIST(I),RMLIST(J),POSION,LATT_A,Q_X,Q_Y,Q_Z,BOND)
           CCOUNT=CCOUNT+1
!           WRITE(IU6,'(I4,2I12,F18.6)') CCOUNT,BALIST(I),RMLIST(J),BOND        
           DBOND(CCOUNT)=BOND
           ATOMR_TMP(CCOUNT,1)=BALIST(I)
           ATOMR_TMP(CCOUNT,2)=RMLIST(J)
         ENDDO
      ENDDO
    
   END SUBROUTINE Bond_ini
  
  SUBROUTINE BSTSELECT(RBOND,ATOMR_TMP,BCOUNT)  
  REAL(q),DIMENSION(7000) ::  RBOND
  INTEGER,DIMENSION(7000,2) :: ATOMR_TMP
  INTEGER :: BCOUNT,I,K
  
  K=0
  DO I=1,NTBS
     IF (AVE_T_BOND(I) .LE. RCUT) THEN
         K=K+1
         RBOND(K)=AVE_T_BOND(I)
         ATOMR_TMP(K,1)=T_RATM(I,1)
         ATOMR_TMP(K,2)=T_RATM(I,2)
         !WRITE(IU6,*)'ATOMR=', ATOMR_TMP(K,:),'TRATM=',T_RATM(I,:)
     ENDIF
  ENDDO
  
  BCOUNT=K
  END SUBROUTINE BSTSELECT

!**********************************************************************
!  BROUTINE FOR OUTPUT INFORMATION in boost steps
!**********************************************************************

    SUBROUTINE BoostOut(NSTEP,NREG,POSION,TIFOR,BOOST_BOND,EPSR_Q,EPSR_MAX,A_EPS_M,SUM_V,BOOST_FACT,&
                     &AVE_BOOST_FACT,SDTIME,SPTIME,SDTIME_B,TADF,ORI_FOR,MI,MATOM1,MATOM2,MFORCE)
                     
      INTEGER :: I,J,k,NSTEP,NBLOCK=200,NREG
      INTEGER :: BI,BATM1,BATM2,MI,MATOM1,MATOM2
      Real(q),DIMENSION(NBBS) :: EPSR_Q,BOOST_BOND
      REAL(q),DIMENSION(3,NIONS) :: POSION,TIFOR,TADF,ORI_FOR
      REAL(q),DIMENSION(3,3) :: LATT_A
      REAL(q) :: EPSR_MAX,A_EPS_M,BOOST_FACT,AVE_BOOST_FACT,SUM_V
      REAL(q) :: DFORCER,MFORCE,VTMP(3)
      REAL(q) :: SDTIME,SPTIME,SDTIME_B
              
      IF (NREG == RMDS+1) THEN
           
           WRITE(IU6,'(/,A,I5,A)')'The Thermal Average of All ',NTBS,' Tagged Bonds involved with Boosted Atoms:'
           WRITE(IU6,'(A)')'-------------------------------------------------------------------'
           WRITE(IU6,'(2X,A,7X,A,7X,A,7X,A)')'Bond','Atom1','Atom2','Bond Length'

              Do I=1,NTBS
                 WRITE(IU6,'(I4,2I12,F18.6)') I,T_RATM(I,1),T_RATM(I,2),AVE_T_BOND(I)
              ENDDO

           WRITE(IU6,'(A,/)')'-------------------------------------------------------------------'
           WRITE(IU6,'(/,A,I4,A)')'  List of All ',NBBS,'  Boosted Bonds in Equalibrium State '
           WRITE(IU6,'(A)')'-------------------------------------------------------------------'
           WRITE(IU6,'(2X,A,7X,A,7X,A,7X,A)')'Bond','Atom1','Atom2','Bond Length'

               Do I=1,NBBS
                  WRITE(IU6,'(I4,2I12,F18.6)') I,ATOMR(I,1),ATOMR(I,2),ORI_BOND(I)
                ENDDO

           WRITE(IU6,'(A,/)')'-------------------------------------------------------------------'       
      ENDIF

          WRITE(IU6,'(/,A)')'Bond will be boosted in this Step.' 
     
      ! IF (MOD(NSTEP,NBLOCK)==0) THEN
      !      WRITE(IU6,'(/,A,I4,A)')'  List of All ',NBBS,'  Boosted Bonds with Length and Epsilon/QRR'
      !      WRITE(IU6,'(A)')'-----------------------------------------------------------------------------------'
      !      WRITE(IU6,'(2X,A,7X,A,7X,A,7X,A,7X,A,7X,A)')'Bond','Atom1','Atom2','Bond Length','Ori length',' Eps/QRR'
      !      DO I=1,NBBS
      !        WRITE(IU6,'(I4,2I12,3F18.6)') I,ATOMR(I,1),ATOMR(I,2),BOOST_BOND(I),ORI_BOND(I),EPSR_Q(I)
      !      ENDDO
      !      WRITE(IU6,'(A,/)')'-----------------------------------------------------------------------------------'
      !  ENDIF

        WRITE(IU6,'(A,F10.6,10X,A,F14.6)') "EPSR_MAX/QRR=",EPSR_MAX,"A[EPSR_MAX]=",A_EPS_M
        IF(EPSR_MAX .LT. 1.0_q) THEN
            WRITE(IU6,'(A,F10.6,10X,A,F10.6,10X,A,F10.6)')"SUM_V=",SUM_V,"BOOST_FACT=",BOOST_FACT,"Temperture=",BTEBEG
            !WRITE(IU6,'(A,F20.10,3X,A,F20.10,3X,A,F20.10)')"Simulation Time=",SDTIME,"Physical Time=",SPTIME,"Simulation Time in Bond-Boost Steps=",SDTIME_B
            WRITE(IU6,'(A,I3,A,F12.8,A,I4,10X,A,I4)') 'Delta_F(',MI,'  )=',MFORCE,'(EPSR_MAX)-------ATOMR1=',MATOM1,'ATOMR2=',MATOM2
            IF (MOD(NSTEP,NBLOCK)==0) THEN
              WRITE(IU6,'(A,/)')"---------------------------------------------------------------------------------------"

              WRITE(IU6,'(A)')'List of Ion Position and Forces'
              WRITE(IU6,'(A)')"***************************************************************************************"

              WRITE(IU6,'(15X,A)') "POSITION                                   CHANGED FORCE"

              DO I=1,NIONS
                 VTMP=POSION(1:3,I)
                 WRITE(IU6,'(3F12.7,10X,3F12.7)') VTMP,(TADF(J,I),J=1,3)
              ENDDO
              WRITE(IU6,'(A)')"***************************************************************************************"

              WRITE(IU6,'(15X,A)') " POSITION                                  ORIGINAL TOTAL FORCE"

              DO I=1,NIONS
                 VTMP=POSION(1:3,I)
                 WRITE(IU6,'(3F12.7,10X,3F12.7)') VTMP,(ORI_FOR(J,I),J=1,3)
              ENDDO
              WRITE(IU6,"(A)") "**************************************************************************************"


              WRITE(IU6,'(15X,A)')" POSITION                                    NEW TOTAL FORCE"
              DO I=1,NIONS
                 VTMP=POSION(1:3,I)
                 WRITE(IU6,'(3F12.7,10X,3F12.7)') VTMP, (TIFOR(J,I),J=1,3)
              END DO
              WRITE(IU6,'(A)') "***************************************************************************************"
            ENDIF
         ELSEIF (EPSR_MAX .GE. 1.0_q) THEN
            DO I=1,NBBS
              IF (EPSR_Q(I) .EQ. EPSR_MAX) THEN
                  BI=I
                  BATM1=ATOMR(I,1)
                  BATM2=ATOMR(I,2)
              ENDIF
            ENDDO
           WRITE(IU6,'(A,I3,A,I4,10X,A,I4)') 'Broken Bond',BI,'-------BATM1=',BATM1,'BATM2=',BATM2
           WRITE(IU6,'(A,F14.4,10X,A,F10.6)')"AVE_BOOST_FACT=",AVE_BOOST_FACT,"Temperture=",BTEBEG
        ENDIF          

    END SUBROUTINE BoostOut

END MODULE bbm
