#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: neb.F,v 1.9 2006-06-26 23:51:07 graeme Exp $
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
!     - rewrite neb code
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
    private
    public :: neb_force,neb_init

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
    SUBROUTINE neb_force(tinfo_nions,posion,toten,force,latt_a,latt_b,iu6)
      INTEGER :: tinfo_nions,iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)

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
      CALLMPI_C( M_sum_d( comm_chain, posion_all(1,1,1), nions*3*images))

      toten_all=0._q
      toten_all(node)=toten
      CALLMPI_C( M_sum_d (comm_chain,toten_all(1),images))

      spring=abs(spring)

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
      call dirkar(nions,dist_prev,dir2car)
      call dirkar(nions,dist_next,dir2car)

      ! apply periodic boundary conditions to these vector differences
      call set_pbc(dist_prev)
      call set_pbc(dist_next)

      ! define the tangent at each image 
      if (node.ne.1) then
         energy_prev_greater = toten_all(node-1).gt.toten_all(node)
      endif
      if (node.ne.images) then
         energy_next_greater = toten_all(node+1).gt.toten_all(node)
      endif
         
      if (node.eq.images.and.energy_last.ne.0) then
         energy_next_greater = energy_last.gt.toten_all(node)
      else if (node.eq.images.and.energy_last.eq.0) then
         ! if the final energy is unknown, assume final energy is lower
         energy_next_greater = .false.
      endif
         
      if (node.eq.1.and.energy_first.ne.0) then
         energy_prev_greater = energy_first.gt.toten_all(node)
      else if (node.eq.1.and.energy_first.eq.0) then
         energy_prev_greater = .false.
      endif
         
      if (iu6>0) write(iu6,*) &
         'NEB: the previous image is higher in energy:', energy_prev_greater
      if (iu6>0) write(iu6,*) &
         'NEB: the next image is higher in energy    :', energy_next_greater
         
      if (energy_prev_greater.neqv.energy_next_greater) then
         ! not at an extrema
         if (energy_prev_greater) then
            if (iu6>0) write(iu6,*) 'NEB: only prev energy greater' 
            tangent = -dist_prev
         else
            if (iu6>0) write(iu6,*) 'NEB: only next energy greater' 
            tangent = dist_next
         endif
      else
         ! at an extrema
         if (iu6>0) write(iu6,*) 'NEB: image is at an extrema' 
         if (node.ne.1) then
            energy_diff_prev = toten_all(node-1)-toten_all(node)
         else
            energy_diff_prev = energy_first-toten_all(node)
         endif
         if (node.ne.images) then
            energy_diff_next=toten_all(node+1)-toten_all(node)
         else
            energy_diff_next=energy_last-toten_all(node)
         endif
         min_diff_energy=min(abs(energy_diff_prev),abs(energy_diff_next))
         max_diff_energy=max(abs(energy_diff_prev),abs(energy_diff_next))
         if (iu6>0) write(iu6,'(A30,2F12.6)') &
            'NEB: diff energy (min, max): ', min_diff_energy,max_diff_energy
         if (energy_diff_prev.gt.energy_diff_next) then
            tangent = dist_next*min_diff_energy - dist_prev*max_diff_energy
         else
            tangent = dist_next*max_diff_energy - dist_prev*min_diff_energy
         endif
      endif

      ! use (old) central difference tangent
      if (ltangentold) then
         tangent = return_unit(dist_next) - return_unit(dist_prev)
      endif

      ! normalize the tangent vector
      call set_unit(tangent)

!======================================================================
! find climbing image
!======================================================================

      if (lclimb) then
         lclimb_local=.true.
         ! use energy of first image if specified in the INCAR
         if(energy_first.ne.0) then
            max_energy_image=0
            max_energy = energy_first
         else
            max_energy_image=1
            max_energy = toten_all(1)
         endif
         do ni=1,images
            if (toten_all(ni).gt.max_energy) then
               max_energy = toten_all(ni)
               max_energy_image = ni
            endif
         enddo
         if(energy_last.ne.0 .and. energy_last.gt.max_energy) then
            max_energy_image = images+1
            max_energy = energy_last
         endif

         if (max_energy_image > 0 .and. max_energy_image < images+1) then
            climbing_image = max_energy_image
         else
            lclimb_local=.false.
            if(iu6>=0) write(iu6) 'NEB: no climbing image found'
         endif
      endif

!======================================================================
! force projections
!======================================================================

      f_perp = force - Vproj(force,tangent)
      f_par = force - f_perp
!      f_spring_par = spring*(Vproj(dist_prev,tangent)+Vproj(dist_next,tangent))
      f_spring_par = spring*(sqrt(sum(dist_next**2))-sqrt(sum(dist_prev**2)))*tangent

!======================================================================
! Double nudging
!   for details, see: Trygubenko and Wales, JCP 120, 2082 (2005)
!   turn on by setting LDNEB=.TRUE.
!
!   The idea here is to add a portion of the perpendicular spring force
!   to stabilize and possibly accelerate convergence of the NEB,
!   particularly with second order optimizers such as BFGS (IBRION 1).
!   This has not been well tested, and users should be warned that this
!   extra force will cause some corner cutting, so that images will not
!   rigorously converge to the minimum energy path.  Convergence of the
!   climbing image, however, will not be affected by double nudging.
!
!   We are using an algorithm which smoothly turns off the double
!   nudging as the NEB converges.  To use the original method,
!   set LDNEBORG=.TRUE.
!
!======================================================================

      if(ldneb) then
         f_spring_vec = spring*(dist_prev+dist_next)
         f_spring_perp = f_spring_vec-Vproj(f_spring_vec,tangent)
         f_dneb = f_spring_perp-Vproj(f_spring_perp,f_perp)

         ! smoothly turn off the dneb force at convergence
         if (.not.ldneborg) then
            if (sum(f_spring_perp**2) > 0._q) &
               f_dneb = f_dneb*2/PI*atan(sum(f_perp**2)/sum(f_spring_perp**2))
         endif
      else
         f_dneb = 0._q
      endif

!======================================================================
! add up the forces
!======================================================================

      if (lclimb_local .and. node==climbing_image) then
         force = force - 2._q*Vproj(force,tangent)
      else
         force = f_spring_par + f_perp + f_dneb
      endif

!======================================================================
! print neb information
!======================================================================

 4696 format( /' NEB: Tangent'/ &
              ' ----------------------------------------------')
 4697 format(/1x,'NEB: forces: par spring, perp real, dneb ',3f12.6)
 4698 format(1x,'NEB: distance to prev, next image, angle between ',3f12.6)
 4699 format(1x,'NEB: projections on to tangent (spring, real) ',2f12.6,/)

      if (iu6>=0) then
         write(iu6,4696)
         write(iu6,'(3F13.5)') (tangent(:,ni),ni=1,nions)
         write(iu6,4697) sqrt(sum(f_spring_par**2)),sqrt(sum(f_perp**2)),&
                         sqrt(sum(f_dneb**2))
         write(iu6,4698) sqrt(sum(dist_prev**2)),sqrt(sum(dist_next**2)),&
            acos(min(max(sum(return_unit(dist_prev)* &
            return_unit(dist_next)),-1._q),1._q))*180._q/PI
         write(iu6,4699) sum(f_spring_par*tangent),sum(f_par*tangent)
      endif

#endif

 9999 continue

    END SUBROUTINE neb_force

!**********************************************************************
!
! Initialize the chain (repeated image mode)
! Read the spring constant and the two endpoint images, which are kept
! fixed during the optimization
!
!**********************************************************************

    SUBROUTINE neb_init (T_INFO, IO)
      USE base
      TYPE (in_struct) :: IO
      TYPE (type_info) :: T_INFO

      ! needed only temporarily for reading
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt):: LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics) :: DYN
      INTEGER IDUM,IERR,N,idir,node
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      integer tiu6, tiu0
      tiu6 = IO%IU6
      tiu0 = IO%IU0
      nions = T_INFO%NIONS

      ! quick return, if we are not running in image mode
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
      if(tiu6>0) write(tiu6,'(A5,A16,F14.6)') 'NEB:','SPRING',spring


      ! read climbing image flag to move maximum energy image to the saddle

      lclimb=.false.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LCLIMB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,lclimb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LCLIMB'' from file INCAR.'
         STOP
      ENDIF
      if(tiu6>0) write(tiu6,'(A5,A16,L7)') 'NEB:','LCLIMB',lclimb

      ! read old tangent flag

      ltangentold=.false.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LTANGENTOLD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ltangentold,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LTANGENTOLD'' from file INCAR.'
         STOP
      ENDIF
      if(tiu6>0) write(tiu6,'(A5,A16,L7)') 'NEB:','LTANGENTOLD',ltangentold

      ! read dneb flag (true = dneb, false = neb)
      
      ldneb=.false.
      ldneborg=.false.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LDNEB'' from file INCAR.'
         STOP
      ENDIF
      if(tiu6>0) write(tiu6,'(A5,A16,L7)') 'NEB:','LDNEB',ldneb

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEBORG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneborg,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
            WRITE(TIU0,*)'Error reading item ''LDNEBORG'' from file INCAR.'
         STOP
      ENDIF
      if(tiu6>0) write(tiu6,'(A5,A16,L7)') 'NEB:','LDNEBORG',ldneborg

      ! read energies of end points

      energy_first=0.0
      energy_last=0.0

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'EFIRST','=','#',';','F', &
     &            IDUM,energy_first,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''EFIRST'' from file INCAR.'
         STOP
      ENDIF
      if(tiu6>0) write(tiu6,'(A5,A16,F14.6)') 'NEB:','EFIRST',energy_first

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'ELAST','=','#',';','F', &
     &            IDUM,energy_last,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''ELAST'' from file INCAR.'
         STOP
      ENDIF
      if(tiu6>0) write(tiu6,'(A5,A16,F14.6)') 'NEB:','ELAST',energy_last

      ! allocate the positions
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
! for moving endpoints
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
      subroutine set_pbc(v1)
      real(q) :: v1(3,nions)
      call kardir(nions,v1,car2dir)
      v1=mod(v1+100.5_q,1._q)-0.5_q
      call dirkar(nions,v1,dir2car)
      end subroutine set_pbc
!======================================================================
! Returns a unit vector along v1
!======================================================================
      function return_unit(V1)
      real(q) :: v1(3,nions)
      real(q),dimension(3,nions) :: return_unit
      return_unit=v1*(1._q/sqrt(sum(v1*v1)))
      end function return_unit
!======================================================================
! Sets V1 to be a unit vector
!======================================================================
      subroutine set_unit(V1)
      real(q) :: v1(3,nions)
      v1=return_unit(v1) 
      end subroutine set_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      function vproj(v1,v2)
      real(q) :: v1(3,nions),v2(3,nions),vproj(3,nions)
      vproj=v2*sum(v1*v2)/sum(v2*v2)
      end function vproj

END MODULE neb
