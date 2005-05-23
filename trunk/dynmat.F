#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: dynmat.F,v 1.2 2005-05-23 19:57:28 graeme Exp $
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

!bpu    REAL(q),ALLOCATABLE,SAVE :: posion_all(:,:,:)
!bpu    REAL(q) :: spring=10

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

      CONTAINS
      SUBROUTINE dynmat_force(nions,posion,toten,force,a,b,iu6)
      
      INTEGER :: nions
      INTEGER :: iu6
      REAL(q) :: posion(3,nions)
      REAL(q) :: force(3,nions),toten
      REAL(q) :: a(3,3),b(3,3)
! local variables
      INTEGER i,j,k,l,flag,current_dof,node,iimages
      INTEGER, SAVE :: ni,nj,total_dof,my_first_dof,my_last_dof,dof_per_node
      REAL(q) :: displacement(3,nions),displacevector(3,1),VTMP(3)
      LOGICAL, SAVE :: first=.TRUE.

      node=1
      iimages=images
      if (iimages.eq.0) iimages=1

      if (IU6>=0) write(iu6,11)
 11   format('DOING DYMAT, READING DISPLACECAR')

#if defined(MPI) || defined(MPI_CHAIN)
      if(images>0) node=comm_chain%node_me
#endif

      if (IU6>=0) write(iu6,12) node,iimages
 12   format('FOR ',i3,' OUT OF ',i3)

!----------------------------------------------------------------------
! write high precision forces and positions for extracting matrix
      if (IU6>=0) then
         WRITE(IU6,172)
         DO J=1,NIONS
            write(IU6,*) 'DYNMAT: Loop at J=',J
            VTMP=POSION(1:3,J)
            write(IU6,*) 'DYNMAT: Before DIRKAR'
            CALL DIRKAR(1,VTMP,A)
            write(IU6,*) 'DYNMAT: After DIRKAR'
            WRITE(IU6,176) VTMP
         ENDDO
      
         WRITE(IU6,272)
         DO J=1,NIONS
            WRITE(IU6,276) (force (I,J),I=1,3)
         ENDDO
      endif

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
      ENDDO
      CLOSE(73)
      
      IF (first) THEN           ! first time in this routine
        first=.FALSE.
         if (IU6>=0) THEN
            WRITE(IU6,13) 
            WRITE(IU6,14) 0,0
            WRITE(IU6,15) toten
         endif
 13      format('DYMAT: ******************************')
 14      format('DYMAT: DISPLACEMENT      ',i3,i3)
 15      format('DYMAT: ENERGY            ',f16.6)
 16      format('DYMAT: ------------------------------')
 17      format('DYMAT: DEGREE OF FREEDOM ',i3,i3,i3)
 18      format('DYMAT: FORCE             ',f16.10)
 19      format('DYMAT: VECTOR            ',3f10.6)

! write some info out
         k=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).ne.0) THEN
                  k=k+1
                  if (IU6>=0) THEN
                     WRITE(IU6,16) 
                     WRITE(IU6,17) i,j,k
                  endif
               ENDIF
            ENDDO
         ENDDO

! count total number of displacements in DISPLACECAR
         total_dof=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).ne.0) THEN
                  total_dof=total_dof+1
               ENDIF
            ENDDO
         ENDDO
! define some variables for doing parallel calculation
         dof_per_node=total_dof/iimages
         my_first_dof=dof_per_node*(node-1)+1
         my_last_dof=my_first_dof+dof_per_node-1
         current_dof=0
         flag=0
         DO i=1,nions
            DO j=1,3
               if (displacement(j,i).ne.0.and.flag==0) THEN
                  current_dof=current_dof+1
                  if (my_first_dof.eq.current_dof) THEN
                     nj=j
                     ni=i
                     flag=1
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ELSE                      ! all other times in this routine

! zero displacement vector
         DO i=1,3
            displacevector(i,1)=0
         ENDDO
! load displacement vector
         displacevector(nj,1)=displacement(nj,ni)
         CALL kardir(1,displacevector,b)
! undo previous displacement
         DO j=1,3
            posion(j,ni)=posion(j,ni)-displacevector(j,1)
         ENDDO
! write some info
         k=0
         DO i=1,nions
            DO j=1,3
               IF (displacement(j,i).ne.0) THEN
                  k=k+1
                  if (IU6>=0) then
                     WRITE(IU6,16) 
                     WRITE(IU6,17) i,j,k
                  endif
               ENDIF
            ENDDO
         ENDDO
! update position in displacement array
         nj=nj+1
         IF (nj.gt.3) THEN
            nj=1
            ni=ni+1
         ENDIF
      ENDIF

! find new displacement         
      DO WHILE (displacement(nj,ni).eq.0.and.ni.le.nions)
         nj=nj+1
         IF (nj.gt.3) THEN
            nj=1
            ni=ni+1
         ENDIF
      ENDDO

      if (ni.gt.nions) goto 100 ! for serial job, we must be done
      
! zero displacement vector
      DO j=1,3
         displacevector(j,1)=0
      ENDDO
! update displacement vector
      displacevector(nj,1)=displacement(nj,ni)
      CALL kardir(1,displacevector,b)
      
! modify positions by new displacement vector
      DO j=1,3
         posion(j,ni)=posion(j,ni)+displacevector(j,1)
      ENDDO

 100  continue
      if (IU6>=0) then
         write(iu6,20) 
      endif
 20   format('DYMAT: LEAVING')

    END SUBROUTINE dynmat_force

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

! needed only temporarily
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt)::       LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics)  :: DYN
      INTEGER     IDUM,IERR,N,idir,node
      CHARACTER*1   CHARAC
      COMPLEX(q)  CDUM  ; LOGICAL  LDUM
!R.S
      integer tiu6, tiu0
        tiu6 = IO%IU6
        tiu0 = IO%IU0

! quick return, if we are not running in image mode
      IF (images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)

      node=COMM_CHAIN%NODE_ME
      CALL MAKE_DIR_APP(node)
      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, &
           IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF

      if(tiu6>=0) write(tiu6,*)'DYNMAT: Setting node=',node

#endif
    END SUBROUTINE dynmat_init

  END MODULE dynmat

