#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: neb.F,v 1.20 2009-02-19 03:44:26 graeme Exp $
!
! Module which implements the elastic band and the nudged
! elastic band method (for references see below)
! module becomes active if IMAGES tag is read from the INCAR files
!
! TODO
!   - movable endpoints
!   - string method?
!   - climbing image convergence criteria (default?)
!   - climb to all saddles (default?)
!   - minimize local minima?
!   - correct treatment of first and last images?
!
!**********************************************************************

  MODULE neb
    USE prec
    USE main_mpi
    USE poscar
    USE lattice

    IMPLICIT NONE
    PRIVATE
    PUBLIC :: neb_step,neb_init

    INTEGER :: nions,iu0,iu6
    REAL(q),ALLOCATABLE :: posion_all(:,:,:),latt_a_all(:,:,:),latt_b_all(:,:,:)
    REAL(q) :: h(3,3)
    REAL(q) :: hinv(3,3)
    REAL(q) :: h_prev(3,3), h_next(3,3), havg_prev(3,3), havg_next(3,3)
    REAL(q) :: hinv_prev(3,3), hinv_next(3,3), havginv_prev(3,3),havginv_next(3,3)
    REAL(q) :: spring,energy_first,energy_last
    REAL(q) :: jacobian
    REAL(q),PARAMETER :: PI = 3.141592653589793238_q
    LOGICAL lclimb,ltangentold,ldneb,ldneborg,cell_flag

!**********************************************************************
!
! Routine to calculate forces between the images on the elastic band
!
!**********************************************************************

  CONTAINS
    SUBROUTINE neb_step(optflag,posion,toten,force,stress,latt_a,latt_b)
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3),stress(3,3)
      LOGICAL :: optflag

      ! local variables
      REAL(q) :: dR_prev(3,nions), dR_next(3,nions)
      REAL(q) :: tangent(3,nions)
      REAL(q) :: f_perp(3,nions),f_par(3,nions),f_spring(3,nions)
      REAL(q) :: f_spring_perp(3,nions),f_spring_par(3,nions)
      REAL(q) :: f_spring_vec(3,nions),f_dneb(3,nions)

      REAL(q) :: toten_all(images)
      INTEGER :: node,ni,i,j
 
      ! variables for tangent
      REAL(q) :: max_diff_energy,min_diff_energy
      REAL(q) :: energy_diff_prev,energy_diff_next
      LOGICAL :: energy_prev_greater,energy_next_greater

      ! variables for climbing image
      REAL(q) :: max_energy
      INTEGER :: climbing_image,max_energy_image
      LOGICAL :: lclimb_local

      ! variables for stress neb
      REAL(q) :: dtot_prev(3,nions+3), dtot_next(3,nions+3)
      REAL(q) :: ftot(3,nions+3), tot_tan(3,nions+3)
      REAL(q) :: ftot_perp(3,nions+3),ftot_par(3,nions+3)
      REAL(q) :: ftot_spring_perp(3,nions+3),ftot_spring_par(3,nions+3)
      REAL(q) :: ftot_spring_vec(3,nions+3),ftot_dneb(3,nions+3)
      REAL(q) :: dcell_prev(3,3), dcell_next(3,3)
      REAL(q) :: strainA(3,3),strainB(3,3)
      REAL(q) :: strain_next(3,3),strain_prev(3,3)
      REAL(q) :: displacement, dR

      optflag = .TRUE.
      IF (images==0) GOTO 9999
!      IF (spring==-1000) RETURN

#if defined(MPI) || defined (MPI_CHAIN)
!
! communicate all positions to all nodes (brute force, but simple)
!
      node=comm_chain%node_me
      !dir2car=latt_a
      !car2dir=latt_b

      posion_all(:,:,1:images)=0
      posion_all(:,:,node)=posion
! moving endpoint version
!      posion_all(:,:,0:images)=0
!      posion_all(:,:,node-1)=posion
! need a moving endpoint version
      CALLMPI_C( M_SUM_d( comm_chain, posion_all(1,1,1), nions*3*images))
      !IF (iu6>0) WRITE(iu6,*) 'posion_all'
      !IF (iu6>0) WRITE(iu6,*) posion_all

      IF (cell_flag) THEN
        ! get all box dim 
        latt_a_all(:,:,1:images)=0._q
        latt_a_all(:,:,node)=latt_a
        latt_b_all(:,:,1:images)=0._q
        latt_b_all(:,:,node)=latt_b
        CALLMPI_C( M_SUM_d( comm_chain, latt_a_all(1,1,1), 3*3*images))
        CALLMPI_C( M_SUM_d( comm_chain, latt_b_all(1,1,1), 3*3*images))
      ENDIF

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

      ! get distace in direct coord. (how much atoms move)
      dR_prev = posion_all(:,:,node-1)-posion_all(:,:,node)
      dR_next = posion_all(:,:,node+1)-posion_all(:,:,node)
      h=latt_a
      hinv=latt_b

      IF (cell_flag) THEN

        ! how much cell changes 
        dcell_prev = (latt_a_all(:,:,node-1) - latt_a_all(:,:,node))
        dcell_next = (latt_a_all(:,:,node+1) - latt_a_all(:,:,node))
        IF (iu6>0) WRITE(iu6,*) 'NEBCELL: dcell_prev:'
        IF (iu6>0) WRITE(iu6,*) dcell_prev
        IF (iu6>0) WRITE(iu6,*) 'NEBCELL: dcell_next:'
        IF (iu6>0) WRITE(iu6,*) dcell_next

        ! average basis 
        h_prev = latt_a_all(:,:,node-1)
        h_next = latt_a_all(:,:,node+1)
        havg_prev=0.5_q*(h_prev+h)
        havg_next=0.5_q*(h_next+h)
        hinv_prev = latt_b_all(:,:,node-1)
        hinv_next = latt_b_all(:,:,node+1)

        IF (iu6>0) THEN
          WRITE(iu6,*) 'NEBCELL: havg_prev:'
          WRITE(iu6,*) havg_prev
          WRITE(iu6,*) 'NEBCELL: havg_next:'
          WRITE(iu6,*) havg_next
        ENDIF

        ! find the inverse of the avegrage cell for conversion
        CALL inv_latt(havg_prev,havginv_prev)
        CALL inv_latt(havg_next,havginv_next)

        ! postions are in direct coordinates, convert to cartesian
        CALL dirkar(nions,dR_prev,havg_prev)
        CALL dirkar(nions,dR_next,havg_next)
        ! apply periodic boundary conditions to these vector differences
        CALL set_pbc(dR_prev,havg_prev,havginv_prev)
        CALL set_pbc(dR_next,havg_next,havginv_next)
      ELSE
        ! postions are in direct coordinates, convert to cartesian
        CALL dirkar(nions,dR_prev,h)
        CALL dirkar(nions,dR_next,h)
        ! apply periodic boundary conditions to these vector differences
        CALL set_pbc(dR_prev,h,hinv)
        CALL set_pbc(dR_next,h,hinv)
      ENDIF

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
         ! if the final energy is unknown, assume final energy is lower
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

      ! calculate the avg strain
      IF (cell_flag) THEN
        strainA = dcell_prev
        CALL sdotB(strainA,hinv_prev)
        strainB = dcell_prev
        CALL sdotB(strainB,hinv)
        strain_prev = 0.5_q*(strainA+strainB)

        strainA = dcell_next
        CALL sdotB(strainA,hinv_next)
        strainB = dcell_next
        CALL sdotB(strainB,hinv)
        strain_next = 0.5_q*(strainA+strainB)
        IF(iu6>0) THEN
          WRITE(iu6,*) "NEBCELL: average strain previous image"
          WRITE(iu6,'(3F13.5)') strain_prev
          WRITE(iu6,*) "NEBCELL: average strain next image"
          WRITE(iu6,'(3F13.5)') strain_next
        ENDIF
      ENDIF

      IF (cell_flag) THEN
        ! remove rotation deg of freedom from cell
        CALL sdotA(stress,h)
        !IF (iu6>=0) THEN
        !  WRITE(iu6,*) "NEBCELL: stress dot A"
        !  WRITE(iu6,'(3F13.5)') stress
        !ENDIF
        stress(2,1)=0._q
        stress(3,1)=0._q
        stress(3,2)=0._q
        ! convert normailized stress back to stress tensor
        CALL sdotB(stress,hinv)
        ! put everything in a big vector 
        DO i=1,3
          DO j=1,3
            dtot_prev(i,j)=jacobian*strain_prev(i,j)
            dtot_next(i,j)=jacobian*strain_next(i,j)
            ftot(i,j)=stress(i,j)/jacobian
          ENDDO
          DO j=1,nions
            dtot_prev(i,j+3)=dR_prev(i,j)
            dtot_next(i,j+3)=dR_next(i,j)
            ftot(i,j+3)=force(i,j)
          ENDDO
        ENDDO
      ENDIF

      IF (cell_flag) THEN
        IF(iu6>0) THEN
          WRITE(iu6,*) "NEBCELL: distance previous image"
          WRITE(iu6,'(3F13.5)') dtot_prev
          WRITE(iu6,*) "NEBCELL: distance next image"
          WRITE(iu6,'(3F13.5)') dtot_next
        ENDIF
      ENDIF

      IF (energy_prev_greater.NEQV.energy_next_greater) THEN
         ! not at an extrema
         IF (energy_prev_greater) THEN
            IF (iu6>0) WRITE(iu6,*) 'NEB: only prev energy greater'
            tangent = -dR_prev
            IF (cell_flag) tot_tan  = -dtot_prev
         ELSE
            IF (iu6>0) WRITE(iu6,*) 'NEB: only next energy greater'
            tangent = dR_next
            IF (cell_flag) tot_tan  = dtot_next
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
            energy_diff_next = toten_all(node+1)-toten_all(node)
         ELSE
            energy_diff_next = energy_last-toten_all(node)
         ENDIF
         min_diff_energy = MIN(ABS(energy_diff_prev),ABS(energy_diff_next))
         max_diff_energy = MAX(ABS(energy_diff_prev),ABS(energy_diff_next))
         IF (iu6>0) WRITE(iu6,'(A30,2F12.6)') &
            'NEB: diff energy (min, max): ', min_diff_energy,max_diff_energy
         IF (energy_diff_prev.GT.energy_diff_next) THEN
            tangent = dR_next*min_diff_energy - dR_prev*max_diff_energy
            IF (cell_flag)  &
              tot_tan = dtot_next*min_diff_energy - dtot_prev*max_diff_energy
         ELSE
            tangent = dR_next*max_diff_energy - dR_prev*min_diff_energy
            IF (cell_flag) &
              tot_tan = dtot_next*max_diff_energy - dtot_prev*min_diff_energy
         ENDIF
      ENDIF
         
      ! use (old) central difference tangent
      IF (ltangentold) THEN
         tangent = return_unit(dR_next,nions) - return_unit(dR_prev,nions)
         IF (cell_flag) &
           tot_tan = return_unit(dtot_next,nions+3) - return_unit(dtot_prev,nions+3)
      ENDIF

      ! normalize the tangent vector
      tangent=return_unit(tangent,nions)
      IF (cell_flag) THEN
        tot_tan=return_unit(tot_tan,nions+3)
        IF (iu6>=0) THEN
          !WRITE(iu6,*) "Regular tangent"
          !WRITE(iu6,'(3F13.5)') tangent
          WRITE(iu6,*) "NEBCELL: Jacobian"
          WRITE(iu6,'(3F13.5)') jacobian
          WRITE(iu6,*) "NEBCELL: tangent"
          WRITE(iu6,'(3F13.5)') tot_tan
        ENDIF
      ENDIF

!======================================================================
! find climbing image
!======================================================================

      IF (lclimb) THEN
         lclimb_local = .TRUE.
         ! use energy of first image if specified in the INCAR
         IF(energy_first.NE.0) THEN
            max_energy_image = 0
            max_energy = energy_first
         ELSE
            max_energy_image = 1
            max_energy = toten_all(1)
         ENDIF
         DO ni=1,images
            IF (toten_all(ni).GT.max_energy) THEN
               max_energy = toten_all(ni)
               max_energy_image = ni
            ENDIF
         ENDDO
         IF (energy_last.NE.0 .AND. energy_last.GT.max_energy) THEN
            max_energy_image = images+1
            max_energy = energy_last
         ENDIF

         IF (max_energy_image > 0 .AND. max_energy_image < images+1) THEN
            climbing_image = max_energy_image
         ELSE
            lclimb_local = .FALSE.
            IF (iu6>=0) WRITE(iu6,*) 'NEB: no climbing image found'
         ENDIF
      ENDIF

!======================================================================
! force projections
!======================================================================

      f_perp = force - vproj(force,tangent,nions)
      f_par = force - f_perp
      f_spring_par = spring*(SQRT(SUM(dR_next**2))-SQRT(SUM(dR_prev**2)))*tangent

      IF (cell_flag) THEN
        ftot_perp = ftot - vproj(ftot,tot_tan,nions+3)
        ftot_par = ftot - ftot_perp
        ftot_spring_par = spring*(SQRT(SUM(dtot_next**2))-SQRT(SUM(dtot_prev**2)))*tot_tan
      ENDIF
!======================================================================
! Double nudging
!   for details, see: Trygubenko and Wales, JCP 120, 2082 (2005)
!   turn on by setting LDNEB = .TRUE.
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
!   Thanks to Craig Plaisance for fixing a bug in vproj
!======================================================================

      IF (ldneb) THEN
         f_spring_vec = spring*(dR_prev+dR_next)
         f_spring_perp = f_spring_vec-vproj(f_spring_vec,tangent,nions)
         f_dneb = f_spring_perp-vproj(f_spring_perp,f_perp,nions)

         IF (cell_flag) THEN
           ftot_spring_vec = spring*(dtot_prev+dtot_next)
           ftot_spring_perp = ftot_spring_vec-vproj(ftot_spring_vec,tot_tan,nions+3)
           ftot_dneb = ftot_spring_perp-vproj(ftot_spring_perp,ftot_perp,nions+3)
         ENDIF

         ! smoothly turn off the dneb force at convergence
         IF (.NOT.ldneborg) THEN
            IF (SUM(f_spring_perp**2) > 0._q) &
               f_dneb = f_dneb*2/PI*atan(SUM(f_perp**2)/SUM(f_spring_perp**2))
            IF (cell_flag) THEN
              IF (SUM(ftot_spring_perp**2) > 0._q) &
                 ftot_dneb = ftot_dneb*2/PI*atan(SUM(ftot_perp**2)/SUM(ftot_spring_perp**2))
            ENDIF
         ENDIF
      ELSE
         f_dneb = 0._q
         IF (cell_flag) ftot_dneb = 0._q
      ENDIF

!======================================================================
! add up the forces
!======================================================================

      IF (lclimb_local .AND. node==climbing_image) THEN
         IF (cell_flag) THEN
           ftot = ftot - 2._q*vproj(ftot,tot_tan,nions+3)
         ELSE
           force = force - 2._q*vproj(force,tangent,nions)
         ENDIF
      ELSE
         force = f_spring_par + f_perp + f_dneb
         IF (cell_flag) ftot = ftot_spring_par + ftot_perp + ftot_dneb
      ENDIF

      IF (cell_flag) THEN

        IF (iu6>=0) THEN
          WRITE(iu6,*) "NEBCELL: force perp"
          WRITE(iu6,'(3F13.5)') ftot_perp
          WRITE(iu6,*) "NEBCELL: springs"
          WRITE(iu6,'(3F13.5)') ftot_spring_par
          WRITE(iu6,*) "NEBCELL: force vector"
          WRITE(iu6,'(3F13.5)') ftot
        ENDIF
      ENDIF



!======================================================================
! put things back into regular vectors 
!======================================================================
      IF (cell_flag) THEN
        DO i=1,3
          DO j=1,3
            stress(i,j)=ftot(i,j)
          ENDDO
          DO j=1,nions
            force(i,j)=ftot(i,j+3)
          ENDDO
        ENDDO
        IF (iu6>=0) THEN
          WRITE(iu6,*) "NEBCELL: stress"
          WRITE(iu6,'(3F13.5)') stress
        ENDIF
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
         WRITE(iu6,4698) SQRT(SUM(dR_prev**2)),SQRT(SUM(dR_next**2)),&
            acos(min(max(SUM(return_unit(dR_prev,nions)* &
            return_unit(dR_next,nions)),-1._q),1._q))*180._q/PI
         WRITE(iu6,4699) SUM(f_spring_par*tangent),SUM(f_par*tangent)
      ENDIF

 4700 format( /' NEBCELL: Cell Tangent'/ &
              ' ----------------------------------------------')
 4701 format(/1x,'NEBCELL: par spring, perp REAL, dneb ',3f12.6)
 4702 format(1x,'NEBCELL: distance to prev, next image, angle between ',3f12.6)
 4705 format(1x,'NEBCELL: projections on tangent (spring, REAL) ',2f12.6,/)

      IF (cell_flag) THEN
        IF (iu6>=0) THEN
           WRITE(iu6,4700)
           WRITE(iu6,'(3F13.5)') tot_tan
           WRITE(iu6,4701) SQRT(SUM(ftot_spring_par**2)),SQRT(SUM(ftot_perp**2)),&
                           SQRT(SUM(ftot_dneb**2))
           WRITE(iu6,4702) SQRT(SUM(dtot_prev**2)),SQRT(SUM(dtot_next**2)),&
              acos(min(max(SUM(return_unit(dtot_prev,nions+3)* &
              return_unit(dtot_next,nions+3)),-1._q),1._q))*180._q/PI
           WRITE(iu6,4705) SUM(ftot_spring_par*tot_tan),SUM(ftot_par*tot_tan)
        ENDIF
      ENDIF


#endif

 9999 CONTINUE
      RETURN
    END SUBROUTINE neb_step

!**********************************************************************
!
! Initialize the chain (repeated image mode)
! Read the spring constant and the two endpoint images, which are kept
! fixed during the optimization
!
!**********************************************************************

    SUBROUTINE neb_init (T_INFO, IO, chain_jacob)
      USE base
      TYPE (in_struct) :: IO
      TYPE (TYPE_info) :: T_INFO
      REAL(q) :: chain_jacob

      ! needed only temporarily for reading
      INTEGER NIOND,NIONPD,NTYPPD,NTYPD,ibrion
      TYPE (latt):: LATT_CUR
      TYPE (TYPE_info) :: T_I
      TYPE (dynamics) :: DYN
      INTEGER IDUM,IERR,N,idir,node
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM
      REAL(q) omega_first,omega_last

      iu6 = IO%IU6
      iu0 = IO%IU0
      nions = T_INFO%NIONS

      ! quick return, if we are not running in image mode
      IF (images==0) RETURN
#if defined(MPI) || defined(MPI_CHAIN)

      ! warn against using IBRION = 2
      ibrion = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IBRION','=','#',';','I', &
     &            ibrion,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*)'Error reading item ''IBRION'' from file INCAR.'
         STOP
      ENDIF
      IF( (ibrion.EQ.2) ) THEN
        IF (IU0>=0) WRITE(IU0,*) 'WARNING: it is not advised to use IBRION=2 with the NEB.'
      ENDIF
      
!     determines of cell should change in NEB
      cell_flag=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LNEBCELL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,cell_flag,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
        IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LNEBCELL'' from file INCAR.'
        STOP
      ENDIF

      spring=-5.
      ! read the spring constant
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'SPRING','=','#',';','F', &
     &            IDUM,spring,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*) 'Error reading item ''SPRING'' from file INCAR.'
         STOP
      ENDIF

      lclimb=.TRUE.
      ! read climbing image flag to move maximum energy image to the saddle
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LCLIMB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,lclimb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LCLIMB'' from file INCAR.'
         STOP
      ENDIF

      ltangentold=.FALSE.
      ! read old tangent flag
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LTANGENTOLD','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ltangentold,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LTANGENTOLD'' from file INCAR.'
         STOP
      ENDIF

      ldneb=.FALSE.
      ldneborg=.FALSE.
      ! read dneb flag (TRUE = dneb, FALSE = neb)
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEB','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneb,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LDNEB'' from file INCAR.'
         STOP
      ENDIF
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'LDNEBORG','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ldneborg,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*)'Error reading item ''LDNEBORG'' from file INCAR.'
         STOP
      ENDIF
      energy_first=0.0
      energy_last=0.0
      ! read energies of end points
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'EFIRST','=','#',';','F', &
     &            IDUM,energy_first,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*)'Error reading item ''EFIRST'' from file INCAR.'
         STOP
      ENDIF
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'ELAST','=','#',';','F', &
     &            IDUM,energy_last,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*)'Error reading item ''ELAST'' from file INCAR.'
         STOP
      ENDIF

      IF (iu6>0) THEN
        WRITE(iu6,'(A5,A16,F14.6)') 'NEB:','SPRING',spring
        WRITE(iu6,'(A5,A16,L7)') 'NEB:','LCLIMB',lclimb
        WRITE(iu6,'(A5,A16,L7)') 'NEB:','LTANGENTOLD',ltangentold
        WRITE(iu6,'(A5,A16,L7)') 'NEB:','LDNEB',ldneb
        WRITE(iu6,'(A5,A16,L7)') 'NEB:','LDNEBORG',ldneborg
        WRITE(iu6,'(A5,A16,L7)') 'NEB:','LNEBCELL',cell_flag
        WRITE(iu6,'(A5,A16,F14.6)') 'NEB:','EFIRST',energy_first
        WRITE(iu6,'(A5,A16,F14.6)') 'NEB:','ELAST',energy_last
      ENDIF

      ! allocate the positions
      ALLOCATE(posion_all(3,t_info%nions,0:images+1))
      ALLOCATE(latt_a_all(3,3,0:images+1))
      ALLOCATE(latt_b_all(3,3,0:images+1))


      ! read 00/POSCAR file, a little bit of fiddling is required
      idir=0
      CALL MAKE_DIR_APP(idir)

      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IU0>=0) WRITE(IU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF
      posion_all(:,:,0)= DYN%POSION
      latt_a_all(:,:,0)= LATT_CUR%A
      latt_b_all(:,:,0)= LATT_CUR%B
      omega_first = LATT_CUR%OMEGA

      ! read images+1/POSCAR file
      idir=images+1
      CALL MAKE_DIR_APP(idir)

      CALL RD_POSCAR_HEAD(LATT_CUR, T_I, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)
      CALL RD_POSCAR(LATT_CUR, T_I, DYN, &
           NIOND,NIONPD, NTYPD,NTYPPD, IO%IU0, IO%IU6)

      IF (T_I%NIONS /= T_INFO%NIONS) THEN
         IF (IU0>=0) WRITE(IU0,*)'ERROR: image mode number of ions wrong'
         STOP
      ENDIF
      posion_all(:,:,images+1)= DYN%POSION
      latt_a_all(:,:,images+1)= LATT_CUR%A
      latt_b_all(:,:,images+1)= LATT_CUR%B
      omega_last = LATT_CUR%OMEGA

      ! default jacobian is avg legnth/atom 
      jacobian=(0.5_q*(omega_first+omega_last)/DBLE(nions))**(1._q/3._q)*sqrt(DBLE(nions))
      ! read the jacobian 
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'JACOBIAN','=','#',';','F', &
     &            IDUM,spring,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IU0>=0) WRITE(IU0,*) 'Error reading item ''JACOBIAN'' from file INCAR.'
         STOP
      ENDIF
      IF (IU0>=0) WRITE(IU0,*) 'Jacobian: ',jacobian
      chain_jacob = jacobian

      node=COMM_CHAIN%NODE_ME
      CALL MAKE_DIR_APP(node)
! for moving endpoints
!      CALL MAKE_DIR_APP(node-1)
#endif

    END SUBROUTINE neb_init

!**********************************************************************
! Vector Functions
!**********************************************************************
!======================================================================
! Recalulates OMEGA and the INVERSE lattice latt_A=d2c and latt_b=c2d
!======================================================================
      SUBROUTINE inv_latt(d2c,c2d)
      REAL(q) :: d2c(3,3),c2d(3,3),vol
      INTEGER :: I,J
      ! cross product from lattice.F
      CALL EXPRO(c2d(1:3,1),d2c(1:3,2),d2c(1:3,3))
      CALL EXPRO(c2d(1:3,2),d2c(1:3,3),d2c(1:3,1))
      CALL EXPRO(c2d(1:3,3),d2c(1:3,1),d2c(1:3,2))
      vol =c2d(1,1)*d2c(1,1)+c2d(2,1)*d2c(2,1) &
     &      +c2d(3,1)*d2c(3,1)

      DO I=1,3
      DO J=1,3
        c2d(I,J)=c2d(I,J)/vol
      ENDDO
      ENDDO
      END SUBROUTINE inv_latt

!======================================================================
! Sets a vector to have the smallest length consistent the the periodic
! boundary conditions.
!======================================================================
      SUBROUTINE set_pbc(v1,d2c,c2d)
      REAL(q) :: v1(3,nions),d2c(3,3),c2d(3,3)
      CALL kardir(nions,v1,c2d)
      v1=mod(v1+100.5_q,1._q)-0.5_q
      CALL dirkar(nions,v1,d2c)
      END SUBROUTINE set_pbc
!======================================================================
! Returns a unit vector along v1
!======================================================================
      FUNCTION return_unit(v1,nmax)
      INTEGER :: nmax
      REAL(q) :: v1(3,nmax)
      REAL(q),dimension(3,nmax) :: return_unit
      IF (SUM(v1*v1) .NE. 0._q) THEN
        return_unit=v1*(1._q/SQRT(SUM(v1*v1)))
      ELSE
        return_unit=v1*0._q
      ENDIF
      END FUNCTION return_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      FUNCTION vproj(v1,v2,nmax)
      INTEGER :: nmax
      REAL(q) :: v1(3,nmax),v2(3,nmax),vproj(3,nmax)
      IF (SUM(v2*v2) .NE. 0._q) THEN
        vproj=v2*SUM(v1*v2)/SUM(v2*v2)
      ELSE
        vproj=v2*0._q
      ENDIF
      END FUNCTION vproj
!======================================================================
! Functions to convert stress tensor
!======================================================================
      SUBROUTINE sdotA(V,BASIS) ! (sigma dot A^T)T
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER I,J,K
      DIMENSION V(3,3),BASIS(3,3),AC(3,3)
      AC=0
      DO J=1,3
      DO I=1,3
         !A(I,J)=AC(I,J)
         DO K=1,3
            AC(I,J)=AC(I,J) + V(I,K)*BASIS(K,J)
         ENDDO
      ENDDO
      ENDDO
      V=AC

      RETURN
      END SUBROUTINE

      SUBROUTINE sdotB(V,BASIS)
      USE prec
      IMPLICIT REAL(q) (A-H,O-Z)
      INTEGER I,J,K
      DIMENSION V(3,3),BASIS(3,3),AC(3,3)
      AC=0
      DO J=1,3
      DO I=1,3
         !A(I,J)=AC(I,J)
         DO K=1,3
            AC(I,J)=AC(I,J) + V(I,K)*BASIS(J,K)
         ENDDO
      ENDDO
      ENDDO
      V=AC

      RETURN
      END SUBROUTINE

END MODULE neb

