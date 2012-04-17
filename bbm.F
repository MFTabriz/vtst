#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: bbm.F,v 1.33 2007/06/06 16:03:09 Liang Zhang Exp $
!
! This module implements the bond boost method.  For more information refer
! to the web page:
!   http://theory.cm.utexas.edu/vtsttools/
! and the articles 
!
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
      INTEGER :: Nions,IU0,IU5,IU6,NBAS,NRMS
      INTEGER :: i,j,n
      INTEGER,ALLOCATABLE,DIMENSION(:) :: BALIST,RMLIST
      REAL(q) :: QRR,PRR,DVMAX,RCUT,RDTM,BPOTIM,BTEBEG

    CONTAINS

!**********************************************************************
! Initialize the bbm
!**********************************************************************

    SUBROUTINE bbm_init(T_INFO,IO)
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO
    
      INTEGER ::  FLAG,ACOUNT  
      ! needed only temporarily for reading
!      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
!      TYPE (latt):: LATT_CUR
!      TYPE (type_info) :: T_I
!      TYPE (dynamics) :: DYN

      IU0 = IO%IU0
      IU5 = IO%IU5
      IU6 = IO%IU6
      Nions = T_info%nions
      
      WRITE(IU6,'(A)') "====================Bond-Boost Method Initializing===================="
      CALL ReadBBMVar()
      
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
     
 !     WRITE(IU6,'(A,I7,A)') "There are",NBAS,"     boosted atoms:"
 !     WRITE(IU6,'(A,<NBAS>I7)') "Boosted Atoms List:",BALIST
 !     WRITE(IU6,'(A,<NRMS>I7)') "Rest Atoms List:",RMLIST
      WRITE(IU6,'(A)') "===============Bond-Boost Method Initialization Finished==============="

    END SUBROUTINE bbm_Init

!**********************************************************************
!  Bond-Boost Method For each ionic step
!**********************************************************************

    SUBROUTINE bbm_Step(OptFlag,POSION,TOTEN,TIFOR,LATT_A,LATT_B)
      LOGICAL :: OptFlag
      INTEGER :: TINFO_NIONS,IO_IU6,I,J,k,ACOUNT,BCOUNT,NBBS,FLAG,NSTEP=1
      INTEGER :: ATOMR_1,ATOMR_2
!      INTEGER,DIMENSION(NBAS) :: RMLIST
      REAL(q),DIMENSION(500) :: ORI_BOND,BOOST_BOND,RBOND
      INTEGER,DIMENSION(500,2) :: ATOMR,ATOMR_Tmp
      Real(q),Allocatable,DIMENSION(:) :: EPSR,EPSR_Q
      REAL(q),DIMENSION(3,NIONS) :: POSION,TIFOR,ADD_FOR,FICH,FICHR
      REAL(q),DIMENSION(3,3) :: LATT_A,LATT_B
      REAL(q) :: TOTEN,BOND,EPSR_MAX,A_EPS_M,BOOST_FACT,SUM_V,RICH_1,RICH_2,Q_X,Q_Y,Q_Z,VTMP(3)
      REAL(q) :: SDTIME=0.0,SPTIME=0.0,SDTIME_B=0.0 
      SAVE SDTIME,SPTIME,SDTIME_B,NSTEP
      
!      WRITE(IU6,*) "*******I AM TESTING******"
!      WRITE(IU6,*) NRMS
!      WRITE(IU6,*) RMLIST(:)
!      WRITE(IU6,*) "*******I AM TESTING******"

            
      WRITE(IU6,'(/,A)')     ' Bond-TESTING: -----------------'
     WRITE(IU6,'(A,I7)') 'I am testing, NO. of ionic step:',NSTEP
     WRITE(IU6,*) 'POTIM=',BPOTIM 
 
      BCOUNT=0
      CALL Bondrecord(NBAS,NRMS,RCUT,BALIST,RMLIST,POSION,LATT_A,IU6,RBOND,ATOMR_TMP,BCOUNT)
      
      IF (NSTEP == 1) THEN
          ORI_BOND=RBOND
      ELSE         
          BOOST_BOND=RBOND
      ENDIF
 
      ATOMR=ATOMR_TMP
      NBBS=BCOUNT        ! number of boosted bond
      WRITE(IU6,*)"NBBS=",NBBS
      ALLOCATE(EPSR(NBBS),EPSR_Q(NBBS))

      EPSR_MAX=0.0
      IF (NSTEP > 1) THEN
        DO I=1,NBBS
           EPSR_Q(I)=(BOOST_BOND(I)-ORI_BOND(I))/(QRR*ORI_BOND(I))

           WRITE(IU6,*) I,"::EPSR_Q=",EPSR_Q(I)       
        ENDDO
        DO I=1,NBBS
           EPSR_MAX=MAX(EPSR_MAX,ABS(EPSR_Q(I)))
        ENDDO
        A_EPS_M=(1.0_q-EPSR_MAX**2)**2/(1.0_q-PRR**2*EPSR_MAX**2)
        SDTIME=SDTIME+BPOTIM
        WRITE(IU6,*) "MAX_EPSR",EPSR_MAX,"AEPS",A_EPS_M        
      ENDIF       
      
      IF(SDTIME .GT. RDTM) THEN
         IF(EPSR_MAX .LT. 1.0_q) THEN
             DO I=1,NBBS
              SUM_V=0.0_q
              SUM_V=SUM_V+DVMAX*(1.0_q-EPSR_Q(I)**2)/REAL(NBBS)
             ENDDO
             BOOST_FACT=A_EPS_M*SUM_V
             WRITE(IU6,*)"SUM_V=",SUM_V,"BOOST_FACT=",BOOST_FACT,"TEBEG=",BTEBEG
             SPTIME=SPTIME+BPOTIM*EXP(BOOST_FACT*1E5_q/BTEBEG/8.63125_q)
             SDTIME_B=SDTIME_B+BPOTIM
             WRITE(iu6,*)"SDTIME=",SDTIME,"SPTIME=",SPTIME,"SDTIME_B=",SDTIME_B
     
     
             DO I=1,NBBS
                WRITE (IU6,'(3I7,F14.6)') I, ATOMR(I,1),ATOMR(I,2),ORI_BOND(I)
             ENDDO


             DO I=1,NIONS
                DO J=1,3
                   FICH(J,I)=0.0_q
                ENDDO
             ENDDO

             DO I=1,NBBS

                IF(EPSR_Q(I) .LT. EPSR_MAX) THEN

                   RICH_1=2.0_q*A_EPS_M*DVMAX*EPSR_Q(I)/QRR/ORI_BOND(I)/REAL(NBBS)

                   WRITE(IU6,*) "RICH-1=",RICH_1

                   ATOMR_1= ATOMR(I,1)
                   ATOMR_2= ATOMR(I,2)
                   CALL BANearest(ATOMR_2,ATOMR_1,POSION,LATT_A,IU6,Q_X,Q_Y,Q_Z,BOND)   
       
                   ADD_FOR(1,I)=Q_X/BOOST_BOND(I)*RICH_1
                   ADD_FOR(2,I)=Q_Y/BOOST_BOND(I)*RICH_1
                   ADD_FOR(3,I)=Q_Z/BOOST_BOND(I)*RICH_1
       
                   FICH(1,ATOMR_1)=FICH(1,ATOMR_1)+ADD_FOR(1,I)
                   FICH(2,ATOMR_1)=FICH(2,ATOMR_1)+ADD_FOR(2,I)
                   FICH(3,ATOMR_1)=FICH(3,ATOMR_1)+ADD_FOR(3,I)
 
                   FICH(1,ATOMR_2)=FICH(1,ATOMR_2)-ADD_FOR(1,I)
                   FICH(2,ATOMR_2)=FICH(2,ATOMR_2)-ADD_FOR(2,I)
                   FICH(3,ATOMR_2)=FICH(3,ATOMR_2)-ADD_FOR(3,I)   
        
                ELSEIF (EPSR_Q(I) .EQ. EPSR_MAX) THEN
                   SUM_V=0.0

                   DO J=1,NBBS
                       SUM_V=SUM_V+DVMAX*(1.0_q-(EPSR_Q(J))**2)/REAL(NBBS)
                   ENDDO

                   WRITE(IU6,*)"SUM_V2=  ",SUM_V

                   RICH_1=2.0_q*A_EPS_M*DVMAX*EPSR_Q(I)/QRR/ORI_BOND(I)/REAL(NBBS)

                   RICH_2=2.0_q*(1.0_q-(EPSR_Q(I))**2)*EPSR_Q(I)*(2.0_q*(1.0_q-PRR**2*(EPSR_Q(I))**2 )&
                         &-PRR**2*(1.0_q-(EPSR_Q(I))**2))**2/ORI_BOND(I)/(1.0_q-PRR**2*(EPSR_Q(I))**2)**2

                   WRITE(IU6,*)"RICH_1=",RICH_1
                   WRITE(IU6,*)"RICH_2=",RICH_2
 
                   CALL BANearest(ATOMR_2,ATOMR_1,POSION,LATT_A,IU6,Q_X,Q_Y,Q_Z,BOND)

                   ADD_FOR(1,I)=Q_X/BOOST_BOND(I)*(RICH_1+SUM_V*RICH_2)
                   ADD_FOR(2,I)=Q_Y/BOOST_BOND(I)*(RICH_1+SUM_V*RICH_2)
                   ADD_FOR(3,I)=Q_Z/BOOST_BOND(I)*(RICH_1+SUM_V*RICH_2)

                   FICH(1,ATOMR_1)=FICH(1,ATOMR_1)+ADD_FOR(1,I)
                   FICH(2,ATOMR_1)=FICH(2,ATOMR_1)+ADD_FOR(2,I)
                   FICH(3,ATOMR_1)=FICH(3,ATOMR_1)+ADD_FOR(3,I)

                   FICH(1,ATOMR_2)=FICH(1,ATOMR_2)-ADD_FOR(1,I)
                   FICH(2,ATOMR_2)=FICH(2,ATOMR_2)-ADD_FOR(2,I)
                   FICH(3,ATOMR_2)=FICH(3,ATOMR_2)-ADD_FOR(3,I)
                ENDIF
            ENDDO

            WRITE(IU6,*)"**************************************************************"

            WRITE(IU6,*) "POSITION                                       ADDED FORCE"

            DO I=1,NIONS
               VTMP=POSION(1:3,I)
               WRITE(IU6,*) VTMP,  (FICH(J,I),J=1,3)
            ENDDO
            WRITE(IU6,*)"***************************************************************"

            WRITE(IU6,*) " POSITION                                  ORIGINAL TOTAL FORCE"

            DO I=1,NIONS
               VTMP=POSION(1:3,I)
               WRITE(IU6,*) VTMP, (TIFOR(J,I),J=1,3)
            ENDDO
            WRITE(IU6,*) "*******************************************************************"

            TIFOR=TIFOR-FICH

            WRITE(IU6,*)" POSITION                                   CHANGED TOTAL FORCE"
            DO I=1,NIONS
              VTMP=POSION(1:3,I)
              WRITE(IU6,*) VTMP, (TIFOR(J,I),J=1,3)
             END DO
           WRITE(IU6,*) "*********************************************************************"

           END IF

       IF(EPSR_MAX .GE. 1.0_q) THEN

           QRR=0.0_q
           WRITE(IU6,*) " *********THE BOOST IS OVER******************** "
           WRITE(IU6,*) "**********************************************"
       END IF
    END IF

    WRITE(IU6,'(/,A)')     ' Bond-TESTING: -----------------'
    NSTEP=NSTEP+1
    END SUBROUTINE bbm_Step

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
     BTEBEG=1.E-4_q
      CALL RDATAB(.true.,'INCAR',IU5,'TEBEG','=','#',';','F', &
     &            IDUM,BTEBEG,CDUM,LDUM,CHARAC,N,1,IERR)
     BPOTIM=0.5_q
      CALL RDATAB(.true.,'INCAR',IU5,'POTIM','=','#',';','F', &
     &            IDUM,BPOTIM,CDUM,LDUM,CHARAC,N,1,IERR)
     QRR =0.8_q
      CALL RDATAB(.true.,'INCAR',IU5,'QRR','=','#',';','F', &
     &            IDUM,QRR,CDUM,LDUM,CHARAC,INint,1,IERR)
     PRR =0.95_q 
      CALL RDATAB(.true.,'INCAR',IU5,'PRR','=','#',';','F', &
     &            IDUM,PRR,CDUM,LDUM,CHARAC,INint,1,IERR)
     RCUT =3.0_q
      CALL RDATAB(.true.,'INCAR',IU5,'RCUT','=','#',';','F', &
     &            IDUM,RCUT,CDUM,LDUM,CHARAC,INint,1,IERR)
     DVMAX =0.3_q 
      CALL RDATAB(.true.,'INCAR',IU5,'DVMAX','=','#',';','F', &
     &            IDUM,DVMAX,CDUM,LDUM,CHARAC,INint,1,IERR)  
     RDTM =300_q 
      CALL RDATAB(.true.,'INCAR',IU5,'RDTM','=','#',';','F', &
     &            IDUM,RDTM,CDUM,LDUM,CHARAC,INint,1,IERR)
      IF (IU6>=0) THEN
        WRITE(IU6,'(/,A)')     ' Bond-Boost: --------------Input Parameters--------------'
        WRITE(IU6,'(A,11X,I7)')    ' Bond-Boost:  NBAS',NBAS
        WRITE(IU6,'(A,9X,<NBAS>I5)')   ' Bond-Boost:  BALIST',BALIST(:)
        WRITE(IU6,'(A,10X,F14.6)') ' Bond-Boost:  POTIM',BPOTIM
        WRITE(IU6,'(A,10X,F14.6)') ' Bond-Boost:  TEBEG',BTEBEG
        WRITE(IU6,'(A,12X,F14.6)') ' Bond-Boost:  PRR',PRR
        WRITE(IU6,'(A,12X,F14.6)') ' Bond-Boost:  QRR',QRR
        WRITE(IU6,'(A,11X,F14.6)') ' Bond-Boost:  RCUT',RCUT
        WRITE(IU6,'(A,9X,F14.6)') ' Bond-Boost:  DV_MAX',DVMAX 
        WRITE(IU6,'(A,F14.6,/)') ' Bond-Boost:  Reguler MD Time',RDTM
      END IF
    
   END SUBROUTINE ReadBBMVar

!=================================================================
!
!Nearest Neighbor between 2 atoms with PBC 
!
!=================================================================
    SUBROUTINE BANearest(SN1,SN2,POSION,L,IU6,Q_X,Q_Y,Q_Z,BOND)
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

!    WRITE(IU6,'(A)')   'TEST Bond-Boost:  BALIST'
    WRITE(IU6,'(A)') '****************I AM TESTING**************'
    WRITE(IU6,'(A,2I7)') 'BOND between',SN1,SN2
!    WRITE(IU6,'(3F14.6)') POSION(:,:)
!    WRITE(IU6,'(3F14.6)') L(:,:)
!    WRITE(IU6,*) D_x,D_y,D_z
    WRITE(IU6,'(A,F14.6)') 'Distance    ',BOND
!    WRITE(IU6,'(A)') '****************I AM TESTING**************'
   
    END SUBROUTINE BANearest

!===================================================================
!Record the Bond length and the label of corresponding bonded atoms
!===================================================================     
    SUBROUTINE Bondrecord(NBAS,NRM,RCUT,BALIST,RMLIST,POSION,LATT_A,IU6,RBOND,ATOMR_TMP,BCOUNT)
   
      INTEGER :: IU6,I,J,BCOUNT,NRM,NBAS,FLAG
      INTEGER,DIMENSION(NBAS) :: RMLIST,BALIST
      REAL(q),DIMENSION(500) ::  RBOND
      INTEGER,DIMENSION(500,2) :: ATOMR_TMP
      REAL(q),DIMENSION(3,NIONS) :: POSION
      REAL(q),DIMENSION(3,3) :: LATT_A
      REAL(q) :: BOND,RCUT,Q_X,Q_Y,Q_Z
 
!Record bonds between Boosted Atoms   
    DO I=1,NBAS
         DO J=I+1,NBAS
             CALL BANearest(BALIST(I),BALIST(J),POSION,LATT_A,IU6,Q_X,Q_Y,Q_Z,BOND)
                  IF(BOND .LT. RCUT) THEN
                   BCOUNT=BCOUNT+1
                   RBOND(BCOUNT)=BOND
                   ATOMR_TMP(BCOUNT,1)=BALIST(I)
                   ATOMR_TMP(BCOUNT,2)=BALIST(J)
                  ENDIF
            ENDDO
      ENDDO

!Record bonds between Boosted Atoms and Rest Atoms
    DO I=1,NBAS
       DO J=1,NRM
           CALL BANearest(BALIST(I),RMLIST(J),POSION,LATT_A,IU6,Q_X,Q_Y,Q_Z,BOND)
                      IF(BOND .LT. RCUT) THEN
                   BCOUNT=BCOUNT+1
                   RBOND(BCOUNT)=BOND
                   ATOMR_TMP(BCOUNT,1)=BALIST(I)
                   ATOMR_TMP(BCOUNT,2)=RMLIST(J)
                      ENDIF
         ENDDO
      ENDDO

   END SUBROUTINE Bondrecord
     
    END MODULE bbm

