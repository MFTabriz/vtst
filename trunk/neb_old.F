#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: neb_old.F,v 1.1 2006-12-29 23:09:32 graeme Exp $
!
! Module which implements the elastic band and the nudged
! elastic band method (for references see below)
! module becomes active if IMAGES tag is read from the INCAR files
!
! Version 1.01, Dec 2005
! Version 1.02, Jan 2006, added double nudging (LDNUDGE)
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
    
    REAL(q),ALLOCATABLE,SAVE :: posion_all(:,:,:)
    REAL(q) :: spring
    REAL(q) :: spring2,efirst,elast,spower
    INTEGER ispring
    LOGICAL lclimb,ltangent,ldneb,ldnebnew

!**********************************************************************
!
! Routine for forces between the images on the elastic band
!
!**********************************************************************

  CONTAINS
    SUBROUTINE neb_force(nions,posion,toten,force,a,b,iu6)
      INTEGER :: nions,iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: a(3,3),b(3,3)
! local variables
      REAL(q) :: x(3),y(3),rmagp,rmagn,distprev(3,nions),distnext(3,nions)
      REAL(q) :: tangent(3,nions),norm,pdotn,rmag
      REAL(q) :: fspdotc,fdotc,frealperp,fsprpar,fsprperp
      REAL(q) :: fspar,fmaxatom,ftemp,ftotmag
      INTEGER node,i,ni,nj

! variables for variable spring constant
      REAL(q) :: sfactor1,sfactor2
      REAL(q) :: utotmax, utotmin,uave,uave1,uave2
      INTEGER im
      REAL(q) :: toten_all(images)
! variables for new tangent
      REAL(q) :: max_energy,min_energy
      REAL(q) :: energy_prev,energy_next
      LOGICAL energy_flag_prev,energy_flag_next
! variables for climbing image
      INTEGER max_node
! variables for double nudging
      REAL(q) :: fperp(3,nions),fspring(3,nions)
      REAL(q) :: fspringperp(3,nions),fspringpar(3,nions),fspringdneb(3,nions)

      IF (images==0) GOTO 9999
      IF (spring==-1000) RETURN

#if defined(MPI) || defined (MPI_CHAIN)
!
! communicate all positions to all nodes (brute force, but simple)
!
      node=comm_chain%node_me

      posion_all(:,:,1:images)=0
      posion_all(:,:,node)=posion
! moving endpoint version
!      posion_all(:,:,0:images)=0
!      posion_all(:,:,node-1)=posion
! need a moving endpoint version
      CALLMPI_C( M_sum_d( comm_chain, posion_all(1,1,1), nions*3*images))

      toten_all=0.0d0
      toten_all(node)=toten
      CALLMPI_C( M_sum_d (comm_chain,toten_all(1),images))

!======================================================================
! Calculate spring constants for variable spring constant NEB
!
! This calculates a spring constant between two images that is 
! dependant on the average energy between the two compared with the
! maximum energy in the system, scaling the spring constant by that
! difference (the higher the energy, the larger the spring constant
! and correspondingly shorter distance between springs)
!======================================================================
      spring=abs(spring)
      spring2=abs(spring2)

      utotmax=-100000
      utotmin=100000
      do im=1,images
         if (toten_all(im).lt.utotmin) utotmin=toten_all(im)
      enddo
      if (efirst.lt.utotmin) utotmin=efirst
      if (elast.lt.utotmin) utotmin=elast
      
      do im = 0,images
         if (im.eq.0) then
            uave=abs(efirst+toten_all(im+1)-2*utotmin)
         elseif (im.eq.images) then
            uave=abs(toten_all(im)+elast-2*utotmin)
         else
            uave=abs(toten_all(im)+toten_all(im+1)-2*utotmin)
         endif
         if (uave.gt.utotmax) utotmax=uave
      enddo

      if (node.ne.1) uave1=abs(toten_all(node)+toten_all(node-1)-2*utotmin)
      if (node.ne.images) uave2=abs(toten_all(node)+toten_all(node+1)-2*utotmin)
      
      if (node.eq.1) uave1=abs(toten_all(node)+efirst-2*utotmin)
      if (node.eq.images) uave2=abs(toten_all(node)+elast-2*utotmin)
      
      if (ispring.eq.1) then
!     spring factors for power scaling
         sfactor2=spring2*(uave2/utotmax)**spower+spring
         sfactor1=spring2*(uave1/utotmax)**spower+spring
      else if (ispring.eq.2) then
!     spring factors for exponential scaling
         sfactor2=spring2*exp((uave2-utotmax)**spower)+spring
         sfactor1=spring2*exp((uave1-utotmax)**spower)+spring
      else if (ispring.eq.4) then
!     spring factors for exponential scaling
         sfactor2=spring2*(exp(uave2-utotmax))**spower+spring
         sfactor1=spring2*(exp(uave1-utotmax))**spower+spring
      else if (ispring.eq.5) then
!     spring factors for power 10 scaling
         sfactor2=spring2*(2**(uave2-utotmax))**spower+spring
         sfactor1=spring2*(2**(uave1-utotmax))**spower+spring
      else if (ispring.eq.3) then
!     spring factors for constant springs
         sfactor1=spring
         sfactor2=spring
      endif

      if (iu6>0) then
         write(iu6,*) 'NEBV Params: ispring,spring,spring2,spower,efirst,elast'
         write(iu6,*) 'NEBV Params: ',ispring,spring,spring2,spower,efirst,elast
      endif

!======================================================================
! Calculate the tangent at each point and the distance to neighboring
! images 
!
! This uses a definition of the tangent defined by the neighboring
! image which is higher in energy.  For a complete description, see
! JCP 113, 9978 (2000)
!======================================================================
      if (iu6>=0) write(iu6,*) 'NEB: DOING NEB!'
      i=node
      norm=0.
      pdotn=0.
      fspdotc=0.
      fdotc=0.
      
!      sfactor1=ABS(spring)
!      sfactor2=ABS(spring)

      do ni=1,nions
         do nj=1,3
! moving endpoint version needed
! need a special rule for the endpoints (just minimize them)
            x(nj)=posion_all(nj,ni,node+1)-posion_all(nj,ni,node)
            ! minimum image convention
            x(nj)=mod(x(nj)+100.5_q,1._q)-0.5_q
            y(nj)=posion_all(nj,ni,node-1)-posion_all(nj,ni,node)
            ! minimum image convention
            y(nj)=mod(y(nj)+100.5_q,1._q)-0.5_q
         enddo
         call  dirkar(1,x,a)
         call  dirkar(1,y,a)
         do nj=1,3
            distnext(nj,ni)=x(nj)
            distprev(nj,ni)=y(nj)
         enddo
      enddo
      rmagp=1./sqrt(sum(distprev**2))
      rmagn=1./sqrt(sum(distnext**2))

      if (ltangent) then

         if (node.ne.1) then
            energy_flag_prev=toten_all(node-1).gt.toten_all(node)
         endif
         if (node.ne.images) then
            energy_flag_next=toten_all(node+1).gt.toten_all(node)
         endif
         
         if (node.eq.images.and.elast.ne.0) then
            energy_flag_next=elast.gt.toten_all(node)
         else if (node.eq.images.and.elast.eq.0) then
            energy_flag_next=.false.
         endif
         
         if (node.eq.1.and.efirst.ne.0) then
            energy_flag_prev=efirst.gt.toten_all(node)
         else if (node.eq.1.and.efirst.eq.0) then
            energy_flag_prev=.false.
         endif
         
         if (iu6>0) then 
            write(iu6,*) 'NEBD: energy_flag_prev,next: ',energy_flag_prev,energy_flag_next
         endif
         
         if (energy_flag_prev.neqv.energy_flag_next) then
            if (energy_flag_prev) then
               if (iu6>0) write(iu6,*) 'NEBD: Only prev higher' 
               tangent=-distprev
            else
               if (iu6>0) write(iu6,*) 'NEBD: Only next higher' 
               tangent=distnext
            endif
         else
            if (iu6>0) write(iu6,*) 'NEBD: Both higher or lower' 
            if (node.ne.1) then
               energy_prev=toten_all(node-1)-toten_all(node)
            else
               energy_prev=efirst-toten_all(node)
            endif
            if (node.ne.images) then
               energy_next=toten_all(node+1)-toten_all(node)
            else
               energy_next=elast-toten_all(node)
            endif
            min_energy=min(abs(energy_prev),abs(energy_next))
            max_energy=max(abs(energy_prev),abs(energy_next))
            if (iu6>0) write(iu6,*) 'NEBD: min, max energy: ',min_energy,max_energy
            if (energy_prev.gt.energy_next) then
               if (iu6>0) write(iu6,*) 'NEBD: Prev gt Next' 
               tangent=distnext*min_energy-distprev*max_energy
            else
               if (iu6>0) write(iu6,*) 'NEBD: Next gt Prev' 
               tangent=distnext*max_energy-distprev*min_energy
            endif
         endif
      else
         tangent=distnext*rmagn-distprev*rmagp
      endif

      norm=1.0/sqrt(sum(tangent**2))
      fdotc=sum(tangent*force)*norm
      pdotn=sum(distnext*distprev)
      fspdotc=-sfactor1*sqrt(sum(distprev**2))+sfactor2*sqrt(sum(distnext**2))

      rmag=rmagn*rmagp
      pdotn=acos(pdotn*rmag)*180.0/3.141592

!======================================================================
90    FORMAT( '  energy of chain is (eV) ',F16.6,' for this image ',F16.6,/ &
           '  tangential force (eV/A) ',F16.6, / &
           '  left and right image ',2F10.6,' A')

100   FORMAT( ' TANGENT     ',35X,'CHAIN-FORCE (eV/Angst)'/ &
           ' ----------------------------------------------', &
           '-------------------------------------')
110   FORMAT((3F13.5,3X,3F14.6))
120   FORMAT( ' ----------------------------------------------', &
           '-------------------------------------')

101   FORMAT( ' Tangent     '/ &
           ' ----------------------------------------------')
102   FORMAT( ' Spring Force     '/ &
           ' ----------------------------------------------')
103   FORMAT( ' Total Projected Force     '/ &
           ' ----------------------------------------------')
111   FORMAT((3F13.5))
121   FORMAT( ' ----------------------------------------------')

      if (iu6>=0) then
         write(iu6,101)
         write(iu6,111) (tangent(:,ni)*norm,ni=1,nions)
      endif

!======================================================================
! find climbing image
!======================================================================
      if (lclimb) then
         max_energy=-1e9
         do im=1,images
            if (toten_all(im).gt.max_energy) then
               max_energy=toten_all(im)
               max_node=im
            endif
         enddo
      endif

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
!   There is also an option which will smoothly turn off the double
!   nudging as the NEB converges (again, this is not well tested)
!   To use this new method, set LDNEBNEW=.TRUE.
!
!======================================================================
      fperp=force-tangent*fdotc*norm
      fspring=sfactor1*distprev+sfactor2*distnext
      fspringpar=sum(fspring*tangent*norm)*tangent*norm
      fspringperp=fspring-fspringpar
      fspringdneb=fspringperp-sum(fspringperp*fperp)*fperp/sum(fperp**2)
!
! Could try limiting the dneb force to the magnitude of fspringpar
! which will prevent corner cutting upon convergence
!
!      if(sum(fspringdneb**2).gt.sum(fspringpar**2)) then
!         fspringdneb=fspringdneb*sqrt(sum(fspringpar**2))/sqrt(sum(fspringdneb**2))
!      endif
!
! This could be better because it turns off smoothly at convergence
      if (ldnebnew) then
         if (sum(fspringperp**2).ge.0.) &
            fspringdneb = fspringdneb*atan(sum(fperp**2)/sum(fspringperp**2))
      endif

!======================================================================
! generate NEB force
!======================================================================
      frealperp=0.
      fsprpar=0.
      fsprperp=0.
      ftotmag=0.
      fmaxatom=0.
      
      do ni=1,nions
         ftemp=0.
         do nj=1,3
            fspar=fspdotc*tangent(nj,ni)*norm
            fsprpar=fsprpar+fspar**2
            frealperp=frealperp+(force(nj,ni)-fdotc*tangent(nj,ni)*norm)**2
            if (lclimb.and.node.eq.max_node) then
               force(nj,ni)=force(nj,ni)-2*fdotc*tangent(nj,ni)*norm
            else
               force(nj,ni)=force(nj,ni)+fspar-fdotc*tangent(nj,ni)*norm
               if (ldneb.or.ldnebnew) then
                  force(nj,ni)=force(nj,ni)+fspringdneb(nj,ni)
               endif
            endif
            ftotmag=ftotmag+force(nj,ni)**2
            ftemp=ftemp+force(nj,ni)**2
         enddo
         if (ftemp.gt.fmaxatom) fmaxatom=ftemp
      enddo

!      write(iu6,*) 'NEB TOTAL FORCE: ',sqrt(ftotmag)
!      if (iu6>=0) then
!         write(iu6,103)
!         write(iu6,111) (force(:,ni),ni=1,nions)
!         write(iu6,121)
!      endif

      fsprpar=sqrt(fsprpar)
      fsprperp=sqrt(fsprperp)
      frealperp=sqrt(frealperp)
      ftotmag=sqrt(ftotmag/real(nions))
      fmaxatom=sqrt(fmaxatom)

      if (iu6>=0) then
         write(iu6,4699) fsprpar,fsprperp, frealperp
         write(iu6,4698) 1./rmagp, 1./rmagn, pdotn
         write(iu6,4697) sfactor1,sfactor2
         write(iu6,4696) fspdotc,fdotc
         write(iu6,4694) fmaxatom,ftotmag
      endif
 4696 format(1x,'NEB: projections on to tangent (spring, real) ',2f12.6)
 4697 format(1x,'NEB: spring constants (prev,next) ',2f12.6)
 4699 format(1x,'NEB: Forces: par spring, perp spring, perp real ',3f12.6)
 4698 format(1x,'NEB: distance to prev, next image, angle between ',3f12.6)
 4694 format(1x,'NEB: Forces: max atom, RMS ',2f12.6)

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

! needed only temporarily
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD
      TYPE (latt):: LATT_CUR
      TYPE (type_info) :: T_I
      TYPE (dynamics) :: DYN
      INTEGER IDUM,IERR,N,idir,node
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM
!R.S
      integer tiu6, tiu0
      tiu6 = IO%IU6
      tiu0 = IO%IU0

! quick return, if we are not running in image mode
      IF (images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)

! read the spring constant

      spring=-5.

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SPRING','=','#',';','F', &
     &            IDUM,spring,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''SPRING'' from file INCAR.'
         STOP
      ENDIF

! read in spring 2 for variable NEB

      spring2=spring

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SPRING2','=','#',';','F', &
     &            IDUM,spring2,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''SPRING2'' from file INCAR.'
         STOP
      ENDIF

! read climbing image flag to move maximum energy image to the saddle

      lclimb=.false.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LCLIMB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,lclimb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''LCLIMB'' from file INCAR.'
         STOP
      ENDIF

! read tangent flag (true = new, false = old)

      ltangent=.true.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LTANGENT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ltangent,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''LTANGENT'' from file INCAR.'
         STOP
      ENDIF

! read dneb flag (true = dneb, false = neb)
      
      ldneb=.false.
      
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''LDNEB'' from file INCAR.'
         STOP
      ENDIF

      ldnebnew=.false.

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEBNEW','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldnebnew,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''LDNEBNEW'' from file INCAR.'
         STOP
      ENDIF

! read energies of end points

      efirst=0.0
      elast=0.0

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'EFIRST','=','#',';','F', &
     &            IDUM,efirst,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''EFIRST'' from file INCAR.'
         STOP
      ENDIF

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'ELAST','=','#',';','F', &
     &            IDUM,elast,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''ELAST'' from file INCAR.'
         STOP
      ENDIF

! read in variable spring constant formula flag and spring power

      ispring=3
      spower=1.0

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'ISPRING','=','#',';','I', &
     &            ispring,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''ISPRING'' from file INCAR.'
         STOP
      ENDIF

      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SPOWER','=','#',';','F', &
     &            IDUM,spower,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. &
           ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''SPOWER'' from file INCAR.'
         STOP
      ENDIF

      if (tiu6>0) then
         write(tiu6,*) 'NEB Params: spring, spring2, efirst, elast, ispring, spower: ',spring,spring2,efirst,elast,ispring,spower
      endif

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

END MODULE neb
