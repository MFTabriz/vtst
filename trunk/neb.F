#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: neb.F,v 1.14 2006-10-11 18:22:16 graeme Exp $
!
! Module which implements the elastic band and the nudged
! elastic band method (for references see below)
! module becomes active if IMAGES tag is read from the INCAR files
!
! Version 1.01, Dec 2005
! Version 1.02, Jan 2006, added double nudging (LDNUDGE)
! Version 1.03, June 2006
! - Removed variable springs
!
!   TODO
!     - optimizer module with cg, qm, lbfgs 
!     - new optimizer flags, reset IBRION, POTIM
!     - reWRITE neb code
!         - movable endpoints
!         - string method?
!         - climbing image convergence criteria (default?)
!         - climb to all saddles (default?)
!         - minimize local minima?
!         - correct treatment of first and last images?
!
!**********************************************************************

  MODULE neb
    USE prec
    USE main_mpi
    USE poscar
    USE lattice

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: neb_force,neb_init

    INTEGER :: nions    
    REAL(q),ALLOCATABLE :: posion_all(:,:,:)
    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: spring,energy_first,energy_last
    REAL(q),PARAMETER :: PI = 3.141592653589793238_q
    LOGICAL lclimb,ltangentold,ldneb,ldneborg

!**********************************************************************
!
! Routine to calculate forces between the images on the elastic band
!
!**********************************************************************

  CONTAINS
    SUBROUTINE neb_force(optflag,tinfo_nions,posion,toten,force,latt_a,latt_b,iu6)
      INTEGER :: tinfo_nions,iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL :: optflag

      ! local variables
      REAL(q) :: dist_prev(3,nions), dist_next(3,nions)
      REAL(q) :: tangent(3,nions)
      REAL(q) :: f_perp(3,nions),f_par(3,nions),f_spring(3,nions)
      REAL(q) :: f_spring_perp(3,nions),f_spring_par(3,nions)
      REAL(q) :: f_spring_vec(3,nions),f_dneb(3,nions)

      REAL(q) :: toten_all(images)
      INTEGER :: node,ni
 
      ! variables for tangent
      REAL(q) :: max_diff_energy,min_diff_energy
      REAL(q) :: energy_diff_prev,energy_diff_next
      LOGICAL :: energy_prev_greater,energy_next_greater

      ! variables for climbing image
      REAL(q) :: max_energy
      INTEGER :: climbing_image,max_energy_image
      LOGICAL :: lclimb_local

      optflag = .TRUE.
      IF (images==0) GOTO 9999
!      IF (spring==-1000) RETURN

#if defined(MPI) || defined (MPI_CHAIN)
!
! communicate all positions to all nodes (brute force, but simple)
!
      node=comm_chain%node_me
      dir2car=latt_a
      car2dir=latt_b

      posion_all(:,:,1:images)=0
      posion_all(:,:,node)=posion
! moving endpoint version
!      posion_all(:,:,0:images)=0
!      posion_all(:,:,node-1)=posion
! need a moving endpoint version
      CALLMPI_C( M_SUM_d( comm_chain, posion_all(1,1,1), nions*3*images))

      toten_all=0._q
      toten_all(node)=toten
      CALLMPI_C( M_SUM_d (comm_chain,toten_all(1),images))

      spring=ABS(spring)

!======================================================================
! Calculate the tangent at each point and the distance to neighboring
! images 
!
! This uses a definition of the tangent defined by the neighboring
! image which is higher in energy.  For a complete description, see
! JCP 113, 9978 (2000)
!======================================================================

      dist_prev = posion_all(:,:,node-1)-posion_all(:,:,node)
      dist_next = posion_all(:,:,node+1)-posion_all(:,:,node)

      ! postions are in direct coordinates, convert to cartesian
      CALL dirkar(nions,dist_prev,dir2car)
      CALL dirkar(nions,dist_next,dir2car)

      ! apply periodic boundary conditions to these vector differences
      CALL set_pbc(dist_prev)
      CALL set_pbc(dist_next)

      ! define the tangent at each image 
      IF (node.NE.1) THEN
         energy_prev_greater = toten_all(node-1).GT.toten_all(node)
      ENDIF
      IF (node.NE.images) THEN
         energy_next_greater = toten_all(node+1).GT.toten_all(node)
      ENDIF
         
      IF (node.EQ.images.AND.energy_last.NE.0) THEN
         energy_next_greater = energy_last.GT.toten_all(node)
      ELSE IF (node.EQ.images.AND.energy_last.EQ.0) THEN
         ! IF the final energy is unknown, asSUMe final energy is lower
         energy_next_greater = .FALSE.
      ENDIF
         
      IF (node.EQ.1.AND.energy_first.NE.0) THEN
         energy_prev_greater = energy_first.GT.toten_all(node)
      ELSE IF (node.EQ.1.AND.energy_first.EQ.0) THEN
         energy_prev_greater = .FALSE.
      ENDIF
         
      IF (iu6>0) WRITE(iu6,*) &
         'NEB: the previous image is higher in energy:', energy_prev_greater
      IF (iu6>0) WRITE(iu6,*) &
         'NEB: the next image is higher in energy    :', energy_next_greater
         
      IF (energy_prev_greater.NEQV.energy_next_greater) THEN
         ! not at an extrema
         IF (energy_prev_greater) THEN
            IF (iu6>0) WRITE(iu6,*) 'NEB: only prev energy greater' 
            tangent = -dist_prev
         ELSE
            IF (iu6>0) WRITE(iu6,*) 'NEB: only next energy greater' 
            tangent = dist_next
         ENDIF
      ELSE
         ! at an extrema
         IF (iu6>0) WRITE(iu6,*) 'NEB: image is at an extrema' 
         IF (node.NE.1) THEN
            energy_diff_prev = toten_all(node-1)-toten_all(node)
         ELSE
            energy_diff_prev = energy_first-toten_all(node)
         ENDIF
         IF (node.NE.images) THEN
            energy_diff_next=toten_all(node+1)-toten_all(node)
         ELSE
            energy_diff_next=energy_last-toten_all(node)
         ENDIF
         min_diff_energy=MIN(ABS(energy_diff_prev),ABS(energy_diff_next))
         max_diff_energy=MAX(ABS(energy_diff_prev),ABS(energy_diff_next))
         IF (iu6>0) WRITE(iu6,'(A30,2F12.6)') &
            'NEB: diff energy (min, max): ', min_diff_energy,max_diff_energy
         IF (energy_diff_prev.GT.energy_diff_next) THEN
            tangent = dist_next*min_diff_energy - dist_prev*max_diff_energy
         ELSE
            tangent = dist_next*max_diff_energy - dist_prev*min_diff_energy
         ENDIF
      ENDIF

      ! USE (old) central difference tangent
      IF (ltangentold) THEN
         tangent = RETURN_unit(dist_next) - RETURN_unit(dist_prev)
      ENDIF

      ! normalize the tangent vector
      CALL set_unit(tangent)

!======================================================================
! find climbing image
!======================================================================

      IF (lclimb) THEN
         lclimb_local=.TRUE.
         ! USE energy of first image IF specIFied in the INCAR
         IF(energy_first.NE.0) THEN
            max_energy_image=0
            max_energy = energy_first
         ELSE
            max_energy_image=1
            max_energy = toten_all(1)
         ENDIF
         DO ni=1,images
            IF (toten_all(ni).GT.max_energy) THEN
               max_energy = toten_all(ni)
               max_energy_image = ni
            ENDIF
         ENDDO
         IF(energy_last.NE.0 .AND. energy_last.GT.max_energy) THEN
            max_energy_image = images+1
            max_energy = energy_last
         ENDIF

         IF (max_energy_image > 0 .AND. max_energy_image < images+1) THEN
            climbing_image = max_energy_image
         ELSE
            lclimb_local=.FALSE.
            IF(iu6>=0) WRITE(iu6) 'NEB: no climbing image found'
         ENDIF
      ENDIF

!======================================================================
! force projections
!======================================================================

      f_perp = force - vproj(force,tangent)
      f_par = force - f_perp
!      f_spring_par = spring*(vproj(dist_prev,tangent)+vproj(dist_next,tangent))
      f_spring_par = spring*(SQRT(SUM(dist_next**2))-SQRT(SUM(dist_prev**2)))*tangent

!======================================================================
! Double nudging
!   for details, see: Trygubenko and Wales, JCP 120, 2082 (2005)
!   turn on by setting LDNEB=.TRUE.
!
!   The idea here is to add a portion of the perpENDicular spring force
!   to stabilize and possibly accelerate convergence of the NEB,
!   particularly with second order optimizers such as BFGS (IBRION 1).
!   This has not been well tested, and USErs should be warned that this
!   extra force will caUSE some corner cutting, so that images will not
!   rigorously converge to the minimum energy path.  Convergence of the
!   climbing image, however, will not be affected by double nudging.
!
!   We are using an algorithm which smoothly turns off the double
!   nudging as the NEB converges.  To USE the original method,
!   set LDNEBORG=.TRUE.
!
!======================================================================

      IF(ldneb) THEN
         f_spring_vec = spring*(dist_prev+dist_next)
         f_spring_perp = f_spring_vec-vproj(f_spring_vec,tangent)
         f_dneb = f_spring_perp-vproj(f_spring_perp,f_perp)

         ! smoothly turn off the dneb force at convergence
         IF (.NOT.ldneborg) THEN
            IF (SUM(f_spring_perp**2) > 0._q) &
               f_dneb = f_dneb*2/PI*atan(SUM(f_perp**2)/SUM(f_spring_perp**2))
         ENDIF
      ELSE
         f_dneb = 0._q
      ENDIF

!======================================================================
! add up the forces
!======================================================================

      IF (lclimb_local .AND. node==climbing_image) THEN
         force = force - 2._q*vproj(force,tangent)
      ELSE
         force = f_spring_par + f_perp + f_dneb
      ENDIF

!======================================================================
! print neb information
!======================================================================

 4696 format( /' NEB: Tangent'/ &
              ' ----------------------------------------------')
 4697 format(/1x,'NEB: forces: par spring, perp REAL, dneb ',3f12.6)
 4698 format(1x,'NEB: distance to prev, next image, angle between ',3f12.6)
 4699 format(1x,'NEB: projections on to tangent (spring, REAL) ',2f12.6,/)

      IF (iu6>=0) THEN
         WRITE(iu6,4696)
         WRITE(iu6,'(3F13.5)') (tangent(:,ni),ni=1,nions)
         WRITE(iu6,4697) SQRT(SUM(f_spring_par**2)),SQRT(SUM(f_perp**2)),&
                         SQRT(SUM(f_dneb**2))
         WRITE(iu6,4698) SQRT(SUM(dist_prev**2)),SQRT(SUM(dist_next**2)),&
            acos(min(max(SUM(RETURN_unit(dist_prev)* &
            RETURN_unit(dist_next)),-1._q),1._q))*180._q/PI
         WRITE(iu6,4699) SUM(f_spring_par*tangent),SUM(f_par*tangent)
      ENDIF

#endif

 9999 CONTINUE
      RETURN   
    END SUBROUTINE neb_force

!**********************************************************************
!
! Initialize the chain (repeated image mode)
! Read the spring constant and the two ENDpoint images, which are kept
! fixed during the optimization
!
!**********************************************************************

    SUBROUTINE neb_init (T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (TYPE_info) :: T_INFO

      ! needed only temporarily for reading
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt):: LATT_CUR
      TYPE (TYPE_info) :: T_I
      TYPE (dynamics) :: DYN
      INTEGER IDUM,IERR,N,idir,node
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      INTEGER tiu6, tiu0
      tiu6 = IO%IU6
      tiu0 = IO%IU0
      nions = T_INFO%NIONS

      ! quick RETURN, IF we are not running in image mode
      IF (images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)

      ! read the spring constant

      spring=-5.

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SPRING','=','#',';','F', &
     &            IDUM,spring,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''SPRING'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,F14.6)') 'NEB:','SPRING',spring


      ! read climbing image flag to move maximum energy image to the saddle

      lclimb=.FALSE.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LCLIMB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,lclimb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LCLIMB'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,L7)') 'NEB:','LCLIMB',lclimb

      ! read old tangent flag

      ltangentold=.FALSE.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LTANGENTOLD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ltangentold,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LTANGENTOLD'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,L7)') 'NEB:','LTANGENTOLD',ltangentold

      ! read dneb flag (TRUE = dneb, FALSE = neb)
      
      ldneb=.FALSE.
      ldneborg=.FALSE.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LDNEB'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,L7)') 'NEB:','LDNEB',ldneb

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEBORG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneborg,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LDNEBORG'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,L7)') 'NEB:','LDNEBORG',ldneborg

      ! read energies of END points

      energy_first=0.0
      energy_last=0.0

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'EFIRST','=','#',';','F', &
     &            IDUM,energy_first,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''EFIRST'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,F14.6)') 'NEB:','EFIRST',energy_first

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'ELAST','=','#',';','F', &
     &            IDUM,energy_last,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''ELAST'' from file INCAR.'
         STOP
      ENDIF
      IF(tiu6>0) WRITE(tiu6,'(A5,A16,F14.6)') 'NEB:','ELAST',energy_last

      ! ALLOCATE the positions
      ALLOCATE(posion_all(3,t_info%nions,0:images+1))

      ! read 00/POSCAR file, a little bit of fiddling is required
      idir=0
      CALL MAKE_DIR_APP(idir)

      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, &
           IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF
      posion_all(:,:,0)= DYN%POSION

      ! read images+1/POSCAR file
      idir=images+1
      CALL MAKE_DIR_APP(idir)

      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, &
           IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IO%IU0>=0) WRITE(TIU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF
      posion_all(:,:,images+1)= DYN%POSION

      node=COMM_CHAIN%NODE_ME
      CALL MAKE_DIR_APP(node)
! for moving ENDpoints
!      CALL MAKE_DIR_APP(node-1)
#endif

    END SUBROUTINE neb_init

!**********************************************************************
!
! Vector Functions
!
!**********************************************************************

!======================================================================
! Sets a vector to have the smallest length consistent the the periodic
! boundary conditions. This should really be changed to use the
! Wigner-Sitz cell
!======================================================================
      SUBROUTINE set_pbc(v1)
      REAL(q) :: v1(3,nions)
      CALL kardir(nions,v1,car2dir)
      v1=mod(v1+100.5_q,1._q)-0.5_q
      CALL dirkar(nions,v1,dir2car)
      END SUBROUTINE set_pbc
!=====e================================================================
! Returns a unit vector along v1
!======================================================================
      FUNCTION return_unit(V1)
      REAL(q) :: v1(3,nions)
      REAL(q),dimension(3,nions) :: return_unit
      return_unit=v1*(1._q/SQRT(SUM(v1*v1)))
      END FUNCTION return_unit
!======================================================================
! Sets V1 to be a unit vector
!======================================================================
      SUBROUTINE set_unit(V1)
      REAL(q) :: v1(3,nions)
      v1=RETURN_unit(v1) 
      END SUBROUTINE set_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      FUNCTION vproj(v1,v2)
      REAL(q) :: v1(3,nions),v2(3,nions),vproj(3,nions)
      vproj=v2*SUM(v1*v2)/SQRT(SUM(v2*v2))
      END FUNCTION vproj

END MODULE neb
