#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: dynmat.F,v 1.9 2007-05-07 06:44:14 graeme Exp $
!
! This module, based on Vasps original chain.F, has been rewritten
! to do dynamical matrix calculations.  It reads in the original
! POSCAR file, as well as DISPLACECAR, which has the displacements
! in the various degrees of freedom needed.
!
!**********************************************************************

  MODULE dynmat
    USE prec
    USE main_mpi
    USE poscar
    USE lattice
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: dynmat_step, dynmat_init

    INTEGER :: nions,iu0,iu6

    CONTAINS

!**********************************************************************
!
!  routine for modifying positions per DISPLACECAR
!  note: no forces are modified, only positions
!
!  This routine modifies the positions for the NEXT force call, so that
!  they reflect the next displacement in the DISPLACECAR file.  The 
!  first time, no displacement is made, and this is the center of the
!  numerical derivative.
!
!**********************************************************************

      SUBROUTINE dynmat_step(optflag,posion,toten,force,a,b)
      
      LOGICAL :: optflag
      REAL(q) :: posion(3,nions)
      REAL(q) :: force(3,nions),toten
      REAL(q) :: a(3,3),b(3,3)
! local variables
      INTEGER i,j,k,l,flag,current_dof,node,iimages
      INTEGER, SAVE :: ni,nj,total_dof,my_first_dof,my_last_dof,dof_per_node
      REAL(q) :: displacement(3,nions),displacevector(3,1),VTMP(3)
      LOGICAL, SAVE :: first=.TRUE.

      ! set optflag to give control to method
      optflag=.FALSE.

      node=1
      iimages=images
      IF (iimages.eq.0) iimages=1

      IF (iu6>=0) WRITE(iu6,11)
 11   FORMAT('DOING DYMAT, READING DISPLACECAR')

#if defined(MPI) || defined(MPI_CHAIN)
      IF (images>0) node=comm_chain%node_me
#endif

      IF (iu6>=0) WRITE(iu6,12) node,iimages
 12   FORMAT('FOR ',i3,' OUT OF ',i3)

!----------------------------------------------------------------------
! write high precision forces and positions for extracting matrix
      IF (iu6>=0) THEN
         WRITE(iu6,172)
         DO J=1,NIONS
            WRITE(iu6,*) 'DYNMAT: Loop at J=',J
            VTMP=POSION(1:3,J)
            WRITE(iu6,*) 'DYNMAT: Before DIRKAR'
            CALL DIRKAR(1,VTMP,A)
            WRITE(iu6,*) 'DYNMAT: After DIRKAR'
            WRITE(iu6,176) VTMP
         END DO
      
         WRITE(iu6,272)
         DO J=1,NIONS
            WRITE(iu6,276) (force (I,J),I=1,3)
         END DO
      END IF

 172  FORMAT( ' HIPREC POSITION (A)    '/ &
     &    ' ----------------------------------------------', &
     &    '-------------------------------------')
 176  FORMAT((3F24.15))
      
 272  FORMAT( ' HIPREC TOTAL-FORCE (eV/A)    '/ &
     &    ' ----------------------------------------------', &
     &    '-------------------------------------')
 276  FORMAT((3F24.15))
!----------------------------------------------------------------------

      OPEN(73,FILE='DISPLACECAR')
      DO I=1,NIONS
         READ(73,*) displacement(1,I),displacement(2,I),displacement(3,I)
      END DO
      CLOSE(73)
      
      IF (first) THEN           ! first time in this routine
        first=.FALSE.
         IF(iu6>=0) THEN
            WRITE(iu6,13) 
            WRITE(iu6,14) 0,0
            WRITE(iu6,15) toten
         END IF
 13      FORMAT('DYMAT: ******************************')
 14      FORMAT('DYMAT: DISPLACEMENT      ',i3,i3)
 15      FORMAT('DYMAT: ENERGY            ',f16.6)
 16      FORMAT('DYMAT: ------------------------------')
 17      FORMAT('DYMAT: DEGREE OF FREEDOM ',i3,i3,i3)
 18      FORMAT('DYMAT: FORCE             ',f16.10)
 19      FORMAT('DYMAT: VECTOR            ',3f10.6)

! WRITE some info out
         k=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).NE.0) THEN
                  k=k+1
                  IF (iu6>=0) THEN
                     WRITE(iu6,16) 
                     WRITE(iu6,17) i,j,k
                  END IF
               END IF
            END DO
         END DO

! count total number of displacements in DISPLACECAR
         total_dof=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).NE.0) THEN
                  total_dof=total_dof+1
               END IF
            END DO
         END DO
! define some variables for doing parallel calculation
         dof_per_node=total_dof/iimages
         my_first_dof=dof_per_node*(node-1)+1
         my_last_dof=my_first_dof+dof_per_node-1
         current_dof=0
         flag=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).NE.0.AND.flag==0) THEN
                  current_dof=current_dof+1
                  IF (my_first_dof.EQ.current_dof) THEN
                     nj=j
                     ni=i
                     flag=1
                  END IF
               END IF
            END DO
         END DO
      ELSE                      ! all other times in this routine

! zero displacement vector
         DO i=1,3
            displacevector(i,1)=0
         END DO
! load displacement vector
         displacevector(nj,1)=displacement(nj,ni)
         CALL kardir(1,displacevector,b)
! undo previous displacement
         DO j=1,3
            posion(j,ni)=posion(j,ni)-displacevector(j,1)
         END DO
! WRITE some info
         k=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).NE.0) THEN
                  k=k+1
                  IF (iu6>=0) THEN
                     WRITE(iu6,16) 
                     WRITE(iu6,17) i,j,k
                  END IF
               END IF
            END DO
         END DO
! update position in displacement array
         nj=nj+1
         IF (nj.GT.3) THEN
            nj=1
            ni=ni+1
         END IF
      END IF

! find new displacement         
      DO WHILE (displacement(nj,ni).EQ.0.AND.ni.LE.nions)
         nj=nj+1
         IF (nj.GT.3) THEN
            nj=1
            ni=ni+1
         END IF
      END DO

      IF (ni.GT.nions) GOTO 100 ! for serial job, we must be done
      
! zero displacement vector
      DO j=1,3
         displacevector(j,1)=0
      END DO
! update displacement vector
      displacevector(nj,1)=displacement(nj,ni)
      CALL kardir(1,displacevector,b)
      
! modify positions by new displacement vector
      DO j=1,3
         posion(j,ni)=posion(j,ni)+displacevector(j,1)
      END DO

 100  CONTINUE

      IF (iu6>=0) WRITE(iu6,20)
 20   FORMAT('DYMAT: LEAVING')

    END SUBROUTINE dynmat_step

!**********************************************************************
!
! initialize the chain (repeated image mode)
! read the spring constant
! and  the two outer images, these images are kept fixed
! during the entire simulation
!
!**********************************************************************

    SUBROUTINE dynmat_init (T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO

      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt) :: LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics)  :: DYN
      INTEGER     IDUM,IERR,N,idir,node
      CHARACTER*1   CHARAC
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM

      nions=T_INFO%nions
      iu6 = IO%IU6
      iu0 = IO%IU0

! quick return, if we are not running in image mode
      IF (images==0) RETURN

#if defined(MPI) || defined(MPI_CHAIN)
      node=COMM_CHAIN%NODE_ME
      CALL MAKE_DIR_APP(node)
      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%iu6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD,IO%IU0, IO%iu6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (iu0>=0) WRITE(iu0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF

      if(iu6>=0) WRITE(iu6,*) 'DYNMAT: Setting node=',node
#endif

    END SUBROUTINE dynmat_init

  END MODULE dynmat

