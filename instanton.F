#include "symbol.inc"
!**********************************************************************!
!
! Version 0.01, May 2006 (AA)
!
!**********************************************************************

  MODULE instanton
    USE base
    USE prec
    USE main_mpi
    USE poscar
    USE lattice

    IMPLICIT NONE

    SAVE
    PRIVATE
    PUBLIC :: instanton_step,instanton_init

    INTEGER nions

#if defined(MPI) || defined(MPI_CHAIN)

! Physical constants
    REAL(q),PARAMETER :: kB=8.61738573e-5_q
    REAL(q),PARAMETER :: hbar=6.46538e-2_q
    REAL(q),PARAMETER :: time_unit=1.018046e-14_q

! Variables private to the module
    TYPE (type_info) :: linfo
    INTEGER :: iu0,iu5,iu6
    INTEGER :: nl,nim,ntot,ndim,ndof,natypes,node,fc=0,itr=0,it=0
    INTEGER,PARAMETER :: iuins=281
    INTEGER,ALLOCATABLE,DIMENSION(:) :: natoms,ifra
    REAL(q) :: dR,ltol,dt,maxmove,tpz,eig,eigold,ftol
    REAL(q),DIMENSION(3,3) :: latt_a,latt_b
    REAL(q),ALLOCATABLE,DIMENSION(:) :: mass,Uim,Usp,Uins,al,bl,dl,el,U0
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: ksp,R,F,premat
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:) :: Fsp,Fim,Rim,w,R0,qq,qqold,F0, &
  &                                         vel,Fold,dold,du,d,Ftmp
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:,:) :: PP
    LOGICAL :: ifcg,new=.TRUE.,converged=.FALSE.,line_step=.TRUE.,        &
  &            do_ins,do_pre

#endif

    CONTAINS

!**********************************************************************
!
! Routine for forces between the images on the elastic band
!
!**********************************************************************

    SUBROUTINE instanton_step(optflag,posion,toten,tifor,ina,inb)
      LOGICAL :: optflag
      REAL(q) :: toten
      REAL(q),DIMENSION(3,nions) :: posion,tifor
      REAL(q),DIMENSION(3,3) :: ina,inb
      REAL(q),EXTERNAL :: rane

#if defined(MPI) || defined(MPI_CHAIN)

! Local variables
      INTEGER :: i,j,im,jm,atom,comp,head,k
      REAL(q),DIMENSION(3,ntot,nim) :: t
      REAL(q),DIMENSION(nim) :: fs0
      REAL(q),DIMENSION(3,ntot) :: dx
      REAL(q) :: s0
      REAL,SAVE :: fmax

      latt_a=ina
      latt_b=inb
      optflag=.FALSE.

! Count force calls
      fc=fc+1

! Here I am ...
      node=comm_chain%node_me

! MPI STUFF :: Before calculating the "force" -> pass stuff to the processors
      CALL pass(toten,posion,tifor)

! Spread the force and energy (from the potential) to the other "imaginary" images
       DO im=1,nim/2
         jm=nim-im+1
         Fim(:,:,jm)=Fim(:,:,im)
         Rim(:,:,jm)=Rim(:,:,im)
         Uim(jm)=Uim(im)
       END DO

! Run the Lanczos algorithm to converge to a min-mode
! -> The original force and coordinates are saved as F0 and R0, respectiavely 
      IF (do_ins) THEN 
        IF (.NOT. converged) THEN
          CALL ins_force()         ! -> add the springs to the potential force
          CALL lanczos()           ! -> converged is set to true in here when the mode is found
        END IF
! If converged, then check to see if the effective force is small enough
! -> This should only be done once (the first time) if CG is used
! -> If QM is being used then line_step is always set as .TRUE.
        IF (converged .AND. line_step) THEN
          CALL effective_force()                                   ! -> Project F0 along w 
          fmax=MAXVAL(ABS(F0))
          IF (iu6 >= 0) THEN
!! Calculate S0 here as well
!            DO im=1,nim
!              head=im+1
!              IF (im == nim) head=1
!! Is Rim the correct one to use here ?!?    -> Yes, it is the most current coordinate.
!              dx=Rim(:,:,head)-Rim(:,:,im)       
!              CALL mic(dx,ntot)
!              fs0(im)=0.0_q
!              k=0
!              DO j=1,ntot                   ! -> Do this in the same order as mass was generated
!                DO i=1,3
!                  IF (linfo%lsfor(i,j)) THEN
!                    k=k+1
!                    fs0(im)=fs0(im)+mass(k)*dx(i,j)**2
!!                    write(iuins,*) 'check loop ',i,j,k,mass(k)
!                  END IF
!                END DO
!              END DO
!            END DO
!            s0=nim*kB*tpz/hbar*sum(fs0)
            s0=zero_mode()
            WRITE(iuins,'(1A11,1ES11.4)') ' check S0 :',s0/hbar
            WRITE(iuins,'(11X,1ES11.4)') 2.0_q*sum(Usp)/(kB*tpz)
            WRITE(iuins,'(1A13)') ' .^.^.^.^.^. '
            WRITE(iuins,'(1A4,1I5,1I8,2X,5F16.10)') ' ut ',itr,fc,SUM(U0),SUM(Usp),fmax,eig,s0
            WRITE(iuins,'(1A13)') ' .^.^.^.^.^. '
            DO im=1,nim/2
              DO i=1,ntot
                IF (ANY(F0(:,i,im) /= 0.0_q)) THEN
                  WRITE(iuins,'(1A4,2I6,3F20.10)') ' F: ',im,i,F0(:,i,im)
                END IF
              END DO
            END DO
! Check for convergence
            IF (fmax < ftol) THEN
              tifor=0.0_q
              converged=.FALSE.                     ! -> Skip over the minimizer and exit.
            END IF
          END IF
        END IF
! When converged, move to a minimizer
        IF (converged) THEN
          IF (ifcg) THEN
            IF (line_step) THEN
              CALL cg()                             ! -> effective_force has already been called
            ELSE                                    ! -> min_step = .TRUE. now
              CALL ins_force()
! This may need to be changed if running prefactor in the same run
! Or maybe not since this part should not be reached if the system has
! converged.
              R0=Rim
              F0=Fim
              U0=Uim
              CALL effective_force()
              CALL cg()
              converged=.FALSE.
              new=.TRUE.
            END IF

!            if (iu6 >= 0) then
!              write(iuins,'(1a9,3f14.8)') ' F0  - N ',F0(:,9,1)
!              write(iuins,'(1a9,3f14.8)') ' F0  - H ',F0(:,10,1)
!              write(iuins,'(a)') ' . '
!              write(iuins,'(1a9,3f14.8)') ' Fold- N ',Fold(:,9,1)
!              write(iuins,'(1a9,3f14.8)') ' Fold- H ',Fold(:,10,1)
!              write(iuins,'(a)') ' . '
!              write(iuins,'(1a9,3f14.8)') ' Ftmp- N ',Ftmp(:,9,1)
!              write(iuins,'(1a9,3f14.8)') ' Ftmp- H ',Ftmp(:,10,1)
!              write(iuins,'(a)') ' . '
!            end if

          ELSE
! -> Run quick-min
          END IF
        END IF
      END IF

      IF (do_pre) THEN
        CALL prefactor(tifor)
      END IF 


! MPI STUFF :: Gather stuff afterwards
!      CALL gather(toten,posion,tifor)
      posion=Rim(:,:,node)
      CALL kardir(ntot,posion,latt_b)

#endif

      RETURN
    END SUBROUTINE instanton_step

!**********************************************************************
!
! Initialize the instanton
!
!**********************************************************************

    SUBROUTINE instanton_init(t_info,io)
      TYPE (in_struct) :: io
      TYPE (type_info) :: t_info
      REAL(q),EXTERNAL :: rane

! Local variables
      INTEGER :: IDUM,IERR,Nint,k1,k2,i,j,im,k
      CHARACTER*1 :: CHARAC
      COMPLEX(q) :: CDUM
      LOGICAL :: LDUM,ertil
      REAL(q) :: RDUM,t

      nions=t_info%nions

#if defined(MPI) || defined(MPI_CHAIN)

! Short cuts
      linfo=t_info
      iu0=io%iu0
      iu5=io%iu5
      iu6=io%iu6

! Quick return, if we are not running in image mode (is this necessary ???)
      IF (images == 0) RETURN
                 
! Open the standard output file for instanton calculations.
      IF (iu6 >= 0) THEN
        OPEN(iuins,FILE=DIR_APP(1:DIR_LEN)//'insout.dat')
        WRITE(iuins,'(A)')  &
  &          '# ut    itr     fc      Pot. En.       Spr. En.         Fmax         eig'
      END IF

! Read in the variables necessary for instanton calculations
! Use here the keyletter "Q" to note variables for instanton calculations 
! ... Eigenvalue tolerance for Lanczos iteration
      ltol=1.0e-2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qltol','=','#',';','F', &
     &            IDUM,ltol,CDUM,LDUM,CHARAC,Nint,1,IERR)

! ... Finite step size, used both for building up the Lanczos matrix and
!     the CG line-search
      dR=1.0e-3_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'QdR','=','#',';','F', &
     &            IDUM,dR,CDUM,LDUM,CHARAC,Nint,1,IERR)

! ... Maximum movement for the system in each step
      maxmove=0.3_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qmaxmove','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,Nint,1,IERR)

! ... Timestep if quick-min is used to minimize
      dt=0.15_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qdt','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,Nint,1,IERR)

! ... Maximum size for the Lanczos matrix
      nl=20
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qnl','=','#',';','I', &
     &            nl,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)

! ... Flag to set the CG method to be used to minimze the instanton.
!     Otherwise quick-min is used.
      ifcg=.TRUE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qifcg','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ifcg,CHARAC,Nint,1,IERR)

! ... Temperature
      tpz=150.0_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qtpz','=','#',';','F', &
     &            IDUM,tpz,CDUM,LDUM,CHARAC,Nint,1,IERR)

! ... Find the instanton
      do_ins=.TRUE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qdo_ins','=','#',';','L', &
     &            IDUM,RDUM,CDUM,do_ins,CHARAC,Nint,1,IERR)

! ... Calculate the prefactor
      do_pre=.FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qdo_pre','=','#',';','L', &
     &            IDUM,RDUM,CDUM,do_pre,CHARAC,Nint,1,IERR)

! ... Force tolerance
      ftol=1.0e-6_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'Qftol','=','#',';','F', &
     &            IDUM,ftol,CDUM,LDUM,CHARAC,Nint,1,IERR)

! Need to allocate the necessary variables. We are using only half of the images
! in the FPI chain so R needs to be made twice as big. 
! ... Set parameters private to this module
      nim=2*images
      ntot=t_info%nions
      ndim=3*ntot
!      ndof=count(t_info%lsfor EQV .TRUE.)
      ndof=0
!GH: fix for pgf90
      DO i = 1,Nions
        DO j = 1,3
          IF (t_info%lsfor(j,i)) ndof=ndof+1
        END DO
      END DO
      natypes=t_info%ntyp
      node=comm_chain%node_me

! Echo input (and other variables)
      IF (iu6 >= 0) THEN
        write(iuins,'(a)') ' ########################################### '
        write(iuins,'(1a11,1i3)')    ' I am node ',node
        write(iuins,'(1a12,1es10.3)') ' ltol    -> ',ltol
        write(iuins,'(1a12,1es10.3)') ' ftol    -> ',ftol
        write(iuins,'(1a12,1es10.3)') ' dR      -> ',dR
        write(iuins,'(1a12,1f9.3)')  ' maxmove -> ',maxmove
        write(iuins,'(1a12,1f9.3)')  ' dt      -> ',dt
        write(iuins,'(1a12,1i4)')    ' nl      -> ',nl
        write(iuins,'(1a12,1l3)')    ' ifcg    -> ',ifcg
        write(iuins,'(1a12,1f9.3)')  ' tpz     -> ',tpz
        write(iuins,'(1a12,1l3)')    ' do_ins  -> ',do_ins
        write(iuins,'(1a12,1l3)')    ' do_pre  -> ',do_pre
        write(iuins,'(1a12,1i4)')    ' nim     -> ',nim
        write(iuins,'(1a12,1i4)')    ' ntot    -> ',ntot
        write(iuins,'(1a12,1i4)')    ' ndim    -> ',ndim
        write(iuins,'(1a12,1i4)')    ' ndof    -> ',ndof
        write(iuins,'(1a12,1i4)')    ' natypes -> ',natypes
        write(iuins,'(a)') ' ########################################### '
      END IF

! ... Allocate variables private to this module
      ALLOCATE(natoms(natypes))
      natoms=t_info%nityp

      ALLOCATE(ksp(3,ntot))
      ksp=0.0_q

      ALLOCATE(Uim(nim),Usp(nim),Uins(nim),U0(nim))
      Uim=0.0_q
      Usp=0.0_q
      Uins=0.0_q
      U0=0.0_q

      ALLOCATE(R(3,ntot),F(3,ntot))
      R=0.0_q
      F=0.0_q

      ALLOCATE(al(nl),bl(nl),dl(nl),el(nl))
      al=0.0_q
      bl=0.0_q
      dl=0.0_q
      el=0.0_q

      ALLOCATE(PP(3,ntot,nim,nl))
      PP=0.0_q

      ALLOCATE(Fsp(3,ntot,nim),Fim(3,ntot,nim),Rim(3,ntot,nim),w(3,ntot,nim),        &
  &            R0(3,ntot,nim),qq(3,ntot,nim),qqold(3,ntot,nim),F0(3,ntot,nim))
      Fsp=0.0_q
      Fim=0.0_q
      Rim=0.0_q
      w=0.0_q
      R0=0.0_q
      qq=0.0_q
      qqold=0.0_q
      F0=0.0_q

! Try to be a bit more clever here what needs to be allocated
!      if (do_ins) then
!      end if
!      IF (do_pre) THEN
        ALLOCATE(premat(ndof*nim,ndof*nim))
        ALLOCATE(ifra(ndof),mass(ndof))
        premat=0.0_q
        ifra=0
        mass=0.0_q
        k=0
! Set the spring constants and mass
!--> Make sure that this is correct for orthongal cells as well as non-orthogonal
        IF (iu6 >= 0) THEN
          WRITE(iuins,'(A)') ' Masses & spring constants for each degree of freedom '
        END IF
        DO j=1,ntot
          DO i=1,3
            IF (t_info%lsfor(i,j)) THEN
              k=k+1
              ifra(k)=3*j-3+i
              mass(k)=t_info%pomass(t_info%ityp(j))
              ksp(i,j)=mass(k)*nim*(kB*tpz/hbar)**2
              IF (iu6 >= 0) THEN
                WRITE(iuins,'(1A6,1I4,2X,1A6,1I2,2X,1A7,1F10.6,2X,1A6,1F10.6)')      &
  &                         ' atom ',j,' comp ',i,' mass= ',mass(k),' ksp= ',ksp(i,j)
              END IF
            END IF
          END DO
        END DO

!      END IF

      IF (ifcg) THEN
        ALLOCATE(Fold(3,ntot,nim),dold(3,ntot,nim),du(3,ntot,nim),d(3,ntot,nim),     &
  &              Ftmp(3,ntot,nim))
        Fold=0.0_q
        dold=0.0_q
        du=0.0_q
        d=0.0_q
        Ftmp=0.0_q
      ELSE
        ALLOCATE(vel(3,ntot,nim))
        vel=0.0_q
      END IF

! Do it like this here (i.e. only on node 1) , otherwise all the images get the same
! random numbers since it done with the same seed on separate processors
      node=comm_chain%node_me
      IF (node == 1) THEN
        INQUIRE(FILE=DIR_APP(1:DIR_LEN)//'MODECAR',EXIST=ertil)
        IF (ertil) THEN
          IF (iu6 >= 0) WRITE(iuins,'(A)') ' Found a MODECAR -> reading it ... '
          OPEN(738,FILE=DIR_APP(1:DIR_LEN)//'MODECAR',STATUS='old',ACTION='read')
          DO im=1,nim
            READ(738,*) (w(:,i,im) , i=1,ntot)
          END DO
          CLOSE(738)  
        ELSE
! If no MODECAR is found, then w is generated randomly
          DO im=1,images
            DO i=1,3
              DO j=1,ntot
                IF (t_info%lsfor(i,j)) THEN
                  w(i,j,im)=rane()-0.5_q
                ELSE
                  w(i,j,im)=0.0_q
                END IF
              END DO
            END DO
            w(:,:,nim-im+1)=w(:,:,im)
          END DO
        END IF
      END IF

#endif

      RETURN
    END SUBROUTINE instanton_init


#if defined(MPI) || defined(MPI_CHAIN)


!**********************************************************************
!
! Apply the Lanczos algorithm to isolate the lowest curvarture mode at
! a given configuration. The current configuration is saved in "R0".
!
!**********************************************************************

    SUBROUTINE lanczos()

! Local variables
      REAL(q),SAVE :: alpha,bta
      INTEGER :: im,jm,i,j

! Variables for LAPACK routine DSYEV
      REAL(q),DIMENSION(nl,nl) :: m
      REAL(q),DIMENSION(3*nl) :: work
      INTEGER :: info

      IF (new) THEN
        new=.FALSE.
        bta=SQRT(SUM(w*w))
        F0=Fim
        R0=Rim
        U0=Uim
        itr=itr+1
        IF (iu6 >= 0) THEN
          WRITE(iuins,'(A,1I3,A)') ' .................... Step ',itr,' .................... '
          WRITE(iuins,'(/,A)') ' Build up the Lanczos matrix '
          WRITE(iuins,'(A)') ' --------------- ' 
        END IF       
        qqold=0.0_q
        qq=w/bta
        it=1
        PP=0.0_q
        PP(:,:,:,it)=qq
        Rim=R0+dR*qq
      ELSE
        w=(Fim-F0)-bta*qqold
        alpha=SUM(qq*w)
        dl(it)=alpha
        w=w-alpha*qq
        bta=SQRT(SUM(w*w))
        el(it)=bta
! Check the eigenvalues
        IF (it > 1) THEN
          m(1:it,1:it)=0.0_q
          DO i=1,it
            m(i,i)=-dl(i)/dR
            IF (i < it) THEN
              m(i,i+1)=-el(i)/dR
              m(i+1,i)=m(i,i+1)
            END IF
          END DO
          al(1:it)=0.0_q
          CALL dsyev('V','U',it,m(1:it,1:it),it,al(1:it),work(1:3*it-1),3*it-1,info)          
          eig=al(1)
          converged=(ABS((eig-eigold)/eigold) < ltol)
! If the maximum size of the matrix has been reached, then exit and print a warning
! to unit=iuins
          IF ((it == nl) .AND. (.NOT. converged)) THEN
            converged=.TRUE.
            IF (iu6 >= 0) THEN
              WRITE(iuins,'(A)') ' WARNING WARNING WARNING '
              WRITE(iuins,'(A)') ' Continuing with a eigenvalue that is not converged '
              WRITE(iuins,'(A)') ' ABS((eig-eigold)/eigold) < ltol '
              WRITE(iuins,'(1F14.8,1A3,1F8.3)')  ABS((eig-eigold)/eigold),' <  ',ltol
              WRITE(iuins,'(A)') ' WARNING WARNING WARNING '
            END IF
          END IF
          IF (iu6 >= 0) THEN
            WRITE(iuins,'(1A6,1I3,1A6,1I3,3F16.8)') ' it   ',it,' info ',info,           &
  &                                                 eig,eigold,abs((eig-eigold)/eigold)
            WRITE(iuins,'(A)') ' --------------- '
          END IF
          eigold=eig
        ELSE
          eigold=-alpha/dR
        END IF
        IF (converged) THEN
! -> Find the mode
          IF (iu6 >= 0) WRITE(iuins,'(/,A)') ' Converged to a min-mode, move to minimizer '

! Test the lowest eigenpair
!          dl=-dl/dR
!          el=-el/dR
!          bl(1)=m(1,1)*dl(1)+m(2,1)*el(1)
!          do i=2,it-1
!            bl(i)=m(i-1,1)*el(i-1)+m(i,1)*dl(i)+m(i+1,1)*el(i)
!          end do
!          bl(it)=m(it-1,1)*el(it-1)+m(it,1)*dl(it)
!          alpha=0.0_q
!          do i=1,it
!            alpha=alpha+m(i,1)*bl(i)
!          end do
!          if (iu6 >= 0) write(iuins,*) alpha,eig,alpha/eig

! Extract the real mode
          w=0.0_q
          DO i=1,it
            w=w+m(i,1)*PP(:,:,:,i)  
          END DO
          w=w/SQRT(SUM(w*w))
          IF (iu6 >= 0) THEN
            OPEN(739,FILE=DIR_APP(1:DIR_LEN)//'NEWMODECAR')
            DO im=1,nim
              WRITE(739,'(3es20.10)') (w(:,i,im) , i=1,ntot)
            END DO
            CLOSE(739)
          END IF
        ELSE
! Do one more lanczos step
          it=it+1
          qqold=qq
          qq=w/bta
          PP(:,:,:,it)=qq
          Rim=R0+dR*qq
        END IF
      END IF
     
      RETURN
    END SUBROUTINE lanczos

!**********************************************************************
!
! Sets up and diagonalizes the prefactor matrix for the instanton.
! As this is set up now then it should be run separately after the
! instanton has been converged in a previous run by setting in the INCAR
! Qdo_ins = .FALSE. & Qdo_pre = .TRUE.
!
!**********************************************************************

    SUBROUTINE prefactor(tifor)
      REAL(q),INTENT(OUT),DIMENSION(3,ntot) :: tifor

! Local variables
      LOGICAL,SAVE :: first=.TRUE.,first_step
      INTEGER,SAVE :: active_atom,active_dof,c_ndof,c_step
      INTEGER :: im,jm,tdof,tatm,shi,shj,i,j,n
      REAL(q) :: t,dtau,dtau2,p1,p2,s0,pi

! Variables for LAPACK routine DSYEV
      REAL(q),DIMENSION(ndof*nim) :: vals
      REAL(q),DIMENSION(3*ndof*nim) :: work
      INTEGER :: info

! Ok, lets start this part !!
! We come in here after calculating the force and energy for the INITIAL point
! --> need to do finite difference around it.

! Imaginary time step size
      dtau=hbar/(kB*tpz*nim)
      dtau2=dtau**2

! If it is the first time in here then go back for a new force
! --> need to identify and which atom and component to move
! --> use ifra
      IF (first) THEN 
        first=.FALSE.
        first_step=.TRUE.
        c_ndof=0                               !-> which ndof are we on 
        c_step=0                               !-> which ionic step are we on
        R0=Rim                                 !-> save the initial configuration
        F0=Fim
        U0=Uim
        IF (iu6 >= 0) THEN
          WRITE(iuins,'(1X,1A4,1X,1I4,1X,1A19,/)') 'Need',2*ndof+1,'steps for prefactor'
          WRITE(iuins,'(1X,1A25)') 'Active degrees of freedom'
          WRITE(iuins,'(3(1A4,2X))') 'dof','atom','comp'
          DO i=1,ndof-1
            WRITE(iuins,'(3(1I4,2X))') ifra(i),CEILING(REAL(ifra(i),q)/3),MOD(ifra(i)-1,3)+1
          END DO
          WRITE(iuins,'(3(1I4,2X),/)') ifra(ndof),CEILING(REAL(ifra(ndof),q)/3),MOD(ifra(ndof)-1,3)+1
        END IF
      END IF

      c_step=c_step+1                                       !--> New step

      IF (first_step) THEN
        first_step=.FALSE.

        c_ndof=c_ndof+1
        IF (c_ndof <= ndof) THEN
          active_dof=MOD(ifra(c_ndof)-1,3)+1
          active_atom=CEILING(REAL(ifra(c_ndof),q)/3)
        END IF
!--> Here we come (except for the first time) in with the force from the "+dR" step,
!--> have what we need to calculate the hessian

        IF (c_ndof > 1) THEN
          i=c_ndof-1
          DO im=1,nim/2
            shi=(im-1)*ndof
            DO j=1,ndof
              tdof=MOD(ifra(j)-1,3)+1
              tatm=CEILING(REAL(ifra(j),q)/3)
              premat(shi+i,shi+j)=-(Fim(tdof,tatm,im)-Ftmp(tdof,tatm,im))/(2.0_q*dR)
            END DO
          END DO
        END IF

! If we have reached the last step, then make the hessian blocks symmetric and scale them,
! else move the system again
        IF (c_ndof > ndof) THEN
          DO im=1,nim/2
            jm=nim-im+1
            shi=ndof*(im-1)
            shj=ndof*(jm-1)                   !-> fill in for the geist-images
            DO i=1,ndof
              premat(shi+i,shi+i)=premat(shi+i,shi+i) * dtau2 / mass(i)
              DO j=i+1,ndof
                premat(shi+i,shi+j)=0.5_q*(premat(shi+i,shi+j)+premat(shi+j,shi+i))
                premat(shi+i,shi+j)=premat(shi+i,shi+j) * dtau2 / (sqrt(mass(i)*mass(j)))
                premat(shi+j,shi+i)=premat(shi+i,shi+j)
              END DO
            END DO
            premat(shj+1:shj+ndof,shj+1:shj+ndof)=premat(shi+1:shi+ndof,shi+1:shi+ndof)
          END DO

! Add contribution from the spring constants
          DO im=1,nim
            shi=ndof*(im-1)
            DO i=1,ndof
              premat(shi+i,shi+i)=premat(shi+i,shi+i)+2.0_q
            END DO
            IF (im == 1) THEN
              DO i=1,ndof
                premat(i,ndof+i)=premat(i,ndof+i)-1.0_q
                premat(i,ndof*(nim-1)+i)=premat(i,ndof*(nim-1)+i)-1.0_q
              END DO
            ELSE IF (im == nim) THEN
              DO i=1,ndof
                premat(shi+i,i)=premat(shi+i,i)-1.0_q
                premat(shi+i,shi-ndof+i)=premat(shi+i,shi-ndof+i)-1.0_q
              END DO
            ELSE
              DO i=1,ndof
                premat(shi+i,shi+ndof+i)=premat(shi+i,shi+ndof+i)-1.0_q
                premat(shi+i,shi-ndof+i)=premat(shi+i,shi-ndof+i)-1.0_q
              END DO
            END IF
          END DO

! Write out the Hessian for the chain
          IF (iu6 >= 0) THEN
            OPEN(511,FILE=DIR_APP(1:DIR_LEN)//'premat.fyl',STATUS='replace',ACTION='write')
            n=ndof*nim
            DO i=1,n
              DO j=1,n
                IF (j == n) THEN
                  WRITE(511,'(1F14.8)',ADVANCE='yes') premat(i,j)
                ELSE
                  WRITE(511,'(1F14.8)',ADVANCE='no') premat(i,j)
                END IF
              END DO
            END DO
            CLOSE(511)
          END IF

! Diagonalize the entire matrix
          IF (iu6 >= 0) WRITE(iuins,*) ' .-.-.-.-. DIAGONALIZING .-.-.-.-. '
          CALL dsyev('V','U',ndof*nim,premat,ndof*nim,vals,work(1:3*ndof*nim-1),3*ndof*nim-1,info)
          IF (iu6 >= 0) WRITE(iuins,*) ' info ',info
          IF (iu6 >= 0) THEN
            n=ndof*nim
! Write out the eigenvalues
            OPEN(511,FILE=DIR_APP(1:DIR_LEN)//'eigvals.fyl',STATUS='replace',ACTION='write')
            WRITE(511,'(1F14.8)') (vals(i) , i=1,n)
            CLOSE(511)
            OPEN(511,FILE=DIR_APP(1:DIR_LEN)//'eigvecs.fyl',STATUS='replace',ACTION='write')
            n=ndof*nim
! Write out the eigenvectors as well ... you never know !!
            DO i=1,n
              DO j=1,n
                IF (j == n) THEN
                  WRITE(511,'(1F14.8)',ADVANCE='yes') premat(i,j)
                ELSE
                  WRITE(511,'(1F14.8)',ADVANCE='no') premat(i,j)
                END IF
              END DO
            END DO
            CLOSE(511)
            WRITE(iuins,'(/,1X,1A24)') 'Five lowest eigenvalues:'
            WRITE(iuins,'(1X,1I3,2X,1F14.8)') (i,vals(i) , i=1,5)
          END IF

! Calculate So here as well (twice the instanton action due to the imaginary time kinetic energy)
          Rim=R0
          s0=zero_mode()

! Need the absolute value for the negative eigenvalue and 
! skip the next one ... it SHOULD be zero but isn't ... lets pretend it is !!!
! If results look weird, then this should be checked.
          vals(1)=-vals(1)
          vals(2)=1.0_q
          pi=4.0_q*ATAN(1.0_q)
          p1=SQRT(s0/(2.0_q*pi*hbar))/(dtau*SQRT(PRODUCT(vals)))
          p2=0.5_q*LOG(s0/(2.0_q*pi*hbar))-LOG(dtau)-0.5_q*SUM(LOG(vals))

          IF (iu6 >= 0) THEN
            WRITE(iuins,'(/,1X,1A25,2X,1F8.3)') 'Temperature             =',tpz
            WRITE(iuins,'(1X,1A25,1X,1ES14.6)') 'So                      =',s0
            WRITE(iuins,'(1X,1A25,1X,1ES14.6)') 'Prefactor [tau^-1]      =',p1
            WRITE(iuins,'(1X,1A25,1X,1ES14.6)') 'Prefactor [s^-1]        =',p1/time_unit
            WRITE(iuins,'(1X,1A25,1X,1ES14.6)') 'Log(Prefactor) [tau^-1] =',p2
            WRITE(iuins,'(1X,1A25,1X,1ES14.6)') 'Log(Prefactor) [s^-1]   =',p2-LOG(time_unit)
          END IF

!          d(1)=-d(1)
!!          d(1)=1.0_q2
!          d(2)=1.0_q2
!! Work with the LN of the prefactor
!          pfac = 0.5_q2 * LOG(S0/(2.0_q2*PI*hbar))                &
!   &           -          LOG(dtau)                               &
!   &           - 0.5_q2 * SUM(LOG(d))                             &
!   &           -         (SUM(Uins)-Uref)/(kB*tpz)                &
!   &           -          LOG(time_unit)
!
!!          WRITE(*,'(1A6,3(2X,F18.12))') ' pfac ',tpz,1000.0_q2/tpz,pfac
!!          write(*,*) ' tmpout ',tpz,1000.0_q2/tpz,log(real(size(d),q2))-0.5_q2*sum(log(d))



! All done here ---> KILL KILL KILL 
          tifor=0.0_q
        ELSE
! Move the system in the active dof and recalculate the force
          IF (iu6 >= 0) THEN
            WRITE(iuins,'(1X,1A11,2X,1I4,2X,1A16,2X,1I4)')                                &
  &                       'Moving atom',active_atom,'& component in -',active_dof 
          END IF
          Rim=R0
          DO im=1,nim
            Rim(active_dof,active_atom,im)=Rim(active_dof,active_atom,im)-dR
          END DO

        END IF

      ELSE
        first_step=.TRUE.
! Here we come in the the force from the "-dR" step ... SAVE IT
        Ftmp=Fim
! Move the system in the active dof and recalculate the force
        IF (iu6 >= 0) THEN
          WRITE(iuins,'(1X,1A11,2X,1I4,2X,1A16,2X,1I4)')                                  &
  &                    'Moving atom',active_atom,'& component in +',active_dof 
          END IF

        Rim=R0
        DO im=1,nim
          Rim(active_dof,active_atom,im)=Rim(active_dof,active_atom,im)+dR
        END DO
      END IF

      RETURN
    END SUBROUTINE prefactor

!**********************************************************************

    FUNCTION zero_mode()
      REAL(q) :: zero_mode

      INTEGER :: im,head,i,j,k
      REAL(q),DIMENSION(3,ntot) :: dx
      REAL(q),DIMENSION(nim) :: fs0

! Calculate S0 here as well
      DO im=1,nim
        head=im+1
        IF (im == nim) head=1
! Is Rim the correct one to use here ?!?    -> Yes, it is the most current coordinate.
        dx=Rim(:,:,head)-Rim(:,:,im)
        CALL mic(dx,ntot)
        fs0(im)=0.0_q
        k=0
        DO j=1,ntot                   ! -> Do this in the same order as mass was generated
          DO i=1,3
            IF (linfo%lsfor(i,j)) THEN
              k=k+1
              fs0(im)=fs0(im)+mass(k)*dx(i,j)**2
!              write(iuins,*) 'check loop ',i,j,k,mass(k),dx(i,j)
            END IF
          END DO
        END DO
      END DO

      zero_mode=nim*kB*tpz/hbar*sum(fs0)

    END FUNCTION zero_mode

!**********************************************************************
!
!**********************************************************************

    SUBROUTINE cg()

! Local variables
      LOGICAL,SAVE :: first=.TRUE.
      REAL(q),PARAMETER :: gamma=0.5_q
      REAL(q) :: a1,a2,b1,cr,step,fr,f2d
      REAL(q),SAVE :: f1d

      IF (first) THEN
        first=.FALSE.
        Fold=0.0_q
        dold=0.0_q
      END IF

      IF (line_step) THEN
! Estimate the curvature along the conjugate direction
        line_step=.FALSE.
        IF (iu6 >= 0) THEN
          WRITE(iuins,'(A)') ' Doing line-min in sub. cg() '
          WRITE(iuins,'(A,1F14.8)')  '   Energy -> ',SUM(U0)
        END IF
! F0 is the force at the local point
! -> OR THE PROJECTED FORCES !!!!
        a1=ABS(SUM(F0*Fold))
        a2=SUM(Fold**2)
        IF (a1 < gamma*a2) THEN
          b1=SUM(F0*(F0-Fold))/a2
        ELSE
          b1=0.0_q
        END IF
        d=F0+b1*dold
        du=d/SQRT(SUM(d**2))
        f1d=SUM(F0*du)
        Ftmp=F0
!        Fold=F0
!        dold=d
        Rim=R0+dR*du
      ELSE
        line_step=.TRUE.
! Minimize along the conjugate direction
        IF (iu6 >= 0) THEN
          WRITE(iuins,'(A)') ' Doing min-step in sub. cg() '
          WRITE(iuins,'(A,1F14.8)')  '   Energy -> ',sum(U0)
        end if
        f2d=SUM(F0*du)
        cr=(f1d-f2d)/dR
        IF (cr < 0.0_q) THEN
          step=maxmove
        ELSE
          fr=0.5_q*(f1d+f2d)
          step=fr/cr
          IF (ABS(step) > maxmove) THEN
            step=SIGN(maxmove,step)
          ELSE
            step=step-0.5_q*dR
          END IF
        END IF
        IF (iu6 >= 0) WRITE(iuins,'(1A15,1F12.6)') ' step length : ',step 
        dold=d
!        Fold=0.5_q*(Ftmp+F0)
        Fold=Ftmp
        Rim=R0+step*du
      END IF

      RETURN
    END SUBROUTINE cg

!**********************************************************************
!
!**********************************************************************

    SUBROUTINE qm()

! Local variables
      INTEGER :: i,j,im
      REAL(q),DIMENSION(3,ntot,nim) :: dv,dRq,step
      REAL(q) :: t,fd,fv

      fd=0.0_q
      fv=0.0_q
      DO im=1,nim
        DO i=1,3
          DO j=1,ntot
            t=vel(i,j,im)*F0(i,j,im)
            IF (t < 0.0_q) THEN
              vel(i,j,im)=0.0_q
            ELSE
              fv=fv+t
            END IF
            fd=fd+F0(i,j,im)**2
          END DO
        END DO
      END DO

      DO im=1,nim
        DO i=1,3
          DO j=1,ntot    
            vel(i,j,im)=F0(i,j,im)*fv/fd
            vel(i,j,im)=vel(i,j,im)+F0(i,j,im)*dt
            step(i,j,im)=vel(i,j,im)*dt
          END DO
        END DO
      END DO
 
      t=SQRT(SUM(step*step))
      IF (t > maxmove) THEN
        step=step/t*maxmove
        DO im=1,nim
          DO i=1,3
            DO j=1,ntot
              Rim(i,j,im)=R0(i,j,im)+step(i,j,im)
            END DO
          END DO
        END DO
      ELSE
        DO im=1,nim
          DO i=1,3
            DO j=1,ntot
              Rim(i,j,im)=R0(i,j,im)+step(i,j,im)
            END DO
          END DO
        END DO
      END IF

      if (iu6 >= 0) then
        write(iuins,'(1a17,1f12.7)') ' Quick-min step: ',sqrt(sum(step**2))
      end if

!      dv=F0*dt
!      where(F0*vel < 0.0_q) vel=0.0_q
!      vel=F0*sum(F0*vel)/sum(F0*F0)
!      vel=vel+dv
!      dRq=vel*dt
!      t=sqrt(sum(dRq*dRq))
!      if (t > maxmove) dRq=dRq/t*maxmove
!      Rim=R0+dRq

      RETURN
    END SUBROUTINE qm

!**********************************************************************
!
!**********************************************************************

    SUBROUTINE effective_force()

      REAL(q),DIMENSION(3,ntot,nim) :: Fpar

      Fpar=w*SUM(w*F0)
      IF (eig < 0.0_q) THEN
        F0=F0-2.0_q*Fpar
      ELSE
        F0=-Fpar
      END IF

      RETURN
    END SUBROUTINE effective_force

!**********************************************************************
!
! Calculates the instanton force, i.e. harmonic springs connects to the
! neighboring images + force from the potentials / number of images. 
!
!**********************************************************************

    SUBROUTINE ins_force()

! Local variables
      INTEGER :: im,jm
      REAL(q) :: tt
      REAL(q),DIMENSION(3,ntot) :: dx,t

!! Spread the force and energy (from the potential) to the other "imaginary" images
!      DO im=1,nim/2
!        jm=nim-im+1
!        Fim(:,:,im)=Fim(:,:,im)/nim
!        Fim(:,:,jm)=Fim(:,:,im)
!        Rim(:,:,jm)=Rim(:,:,im)
!        Uim(jm)=Uim(im)
!      END DO
!      Uim=Uim/nim

      Fim=Fim/nim
      Uim=Uim/nim

! Add force and energy from the temperature dependent springs
! ... Here everything should be in cartesian coord's --> double check that
      Fsp=0.0_q
      Usp=0.0_q
      DO im=1,nim/2-1
        dx=Rim(:,:,im+1)-Rim(:,:,im)
        CALL mic(dx,ntot)
        t=ksp*dx
        Fsp(:,:,im)   = Fsp(:,:,im)   + t
        Fsp(:,:,im+1) = Fsp(:,:,im+1) - t
        tt=SUM(dx*t)
        Usp(im)   = Usp(im)   + tt * 0.25_q
        Usp(im+1) = Usp(im+1) + tt * 0.25_q
      END DO
      DO im=1,nim/2
        jm=nim-im+1
        Fsp(:,:,jm)=Fsp(:,:,im)
        Usp(jm)=Usp(im)
      END DO

!      Uim=Uim*nim
!      Fim=Fim*nim

      Uim=Uim+Usp
      Fim=Fim+Fsp

!      Fim=Fsp

      RETURN
    END SUBROUTINE ins_force

!**********************************************************************

    SUBROUTINE pass(toten,posion,tifor)
      REAL(q) :: toten
      REAL(q),DIMENSION(3,ntot) :: posion,tifor     ! -> nions = ntot

! Local variables
      LOGICAL,SAVE :: first=.TRUE.
      
      Uim=0.0_q
      Uim(node)=toten
! .. pass the energy
      CALLMPI_C(M_sum_d(comm_chain,Uim(1),nim))
      R=posion
      CALL dirkar(ntot,R,latt_a)
      Rim(:,:,1:nim)=0.0_q
      Rim(:,:,node)=R
! .. pass the coordinates
      CALLMPI_C(M_sum_d(comm_chain,Rim(1,1,1),3*ntot*nim))
      F=tifor                       !--> Note that forces are always given in cartesian coord's.
      Fim(:,:,1:nim)=0.0_q
      Fim(:,:,node)=F
! .. pass the force (potential)
      CALLMPI_C(M_sum_d(comm_chain,Fim(1,1,1),3*ntot*nim))
! .. pass w
      IF (first) THEN
        first=.FALSE.
        CALLMPI_C(M_sum_d(comm_chain,w(1,1,1),3*ntot*nim))
      END IF

      RETURN
    END SUBROUTINE pass

!**********************************************************************
!
!    SUBROUTINE gather(toten,posion,tifor)
!      REAL(q) :: toten
!      REAL(q),DIMENSION(3,ntot) :: posion,tifor     ! -> nions = ntot
! 
!      toten=Uim(node)
!      tifor=Fim(:,:,node)
!      posion=Rim(:,:,node)
!      CALL kardir(ntot,posion,latt_b)
!
!      RETURN
!    END SUBROUTINE gather
!
!**********************************************************************
!
!   Calculates the minimum-image-convention. 
!
!**********************************************************************

    SUBROUTINE mic(v,n)
      INTEGER :: n
      REAL(q),INTENT(INOUT),DIMENSION(3,n) :: v

      CALL KARDIR(n,v,latt_b)
      v=MOD(v+100.5_q,1.0_q)-0.5_q
      CALL DIRKAR(n,v,latt_a)

      RETURN
    END SUBROUTINE mic 

!**********************************************************************
!
!   Maintains periodic boundary conditions
!
!**********************************************************************

    SUBROUTINE pbc()


      RETURN
    END SUBROUTINE

!**********************************************************************

    SUBROUTINE test_qm()

! -> Test the quick-min part
      CALL ins_force()
      R0=Rim
      F0=Fim
      CALL qm()

      RETURN
    END SUBROUTINE test_qm

!**********************************************************************

    SUBROUTINE test_cg()

! -> Test the conjugate gradients part
      CALL ins_force()
      R0=Rim
      F0=Fim
      eig=100.0_q
      CALL cg()
      if (line_step) then
        if (iu6 >= 0) then
          write(iuins,*) ' '
          write(iuins,*) 'U1,U2,U2',Uim(1:3)
          write(iuins,'(9x,1a15,1f14.8)')            ' Pot. energy : ',sum(Uim)/nim
          write(iuins,'(7x,1a17,1f14.8)')          ' Spring energy : ',sum(Usp)
          write(iuins,'(1i4,3x,1a17,1f14.8)')   itr,' ... Max Feff : ',maxval(abs(F0))
          write(iuins,'(1i4,1a23,1i8)')  itr,' ... Tot force calls : ',fc
          write(iuins,*) ' '
        end if
      end if

      RETURN
    END SUBROUTINE test_cg

!**********************************************************************

#endif

  END MODULE instanton
