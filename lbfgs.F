#include "symbol.inc"
!**********************************************************************
!
! Module for limited-memory bfgs
!
! - Included is the global treatment for use with the NEB (LGLOBAL) 
! - Included is an option to use a line minimzer (LLINEOPT)
! - Updated to include AutoScale INVCURV by Sam Chill
!
!**********************************************************************

  MODULE lbfgs
    USE prec
    USE main_mpi
    USE lattice

    IMPLICIT NONE
    SAVE
    PRIVATE 
    PUBLIC :: lbfgs_step, lbfgs_init  !CALL opt_init from chain_init

    INTEGER :: nions, iu6, memory, itr, images_local
    REAL(q),ALLOCATABLE ::R(:,:,:), Rold(:,:,:), F(:,:,:), Fold(:,:,:)
    REAL(q),ALLOCATABLE :: rho(:), alpha(:)
    REAL(q),ALLOCATABLE :: direction(:,:,:)
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:,:) :: change_in_G, change_in_R

    REAL(q) :: dir2car(3,3), car2dir(3,3)
    REAL(q) :: finite_step, maxmove, invcurv, damping
    LOGICAL :: fdstep, reset_flag, global, lineopt, autoscale

!**********************************************************************
!
! Limited-memory bfgs method
!
!**********************************************************************
  CONTAINS
    SUBROUTINE lbfgs_step(optflag, posion, toten, force, latt_a, latt_b)

      REAL(q) :: posion(3,nions), toten, force(3,nions)
      REAL(q) :: latt_a(3,3), latt_b(3,3)
      LOGICAL optflag, maxmove_flag

      REAL(q) :: beta
      REAL(q) :: a1, a2, curvature, direction_norm
      REAL(q) :: fp1, fp2, Favg, step_size
      INTEGER :: node, bound, im, jm

      maxmove_flag = .FALSE.

      IF(iu6>0) WRITE(iu6,*) "OPT: LBFGS start"
!#if defined(MPI) || defined (MPI_CHAIN)
      ! Convert the position into Cartesian coordinates
      dir2car = latt_a
      car2dir = latt_b
      CALL dirkar(nions, posion, latt_a)
      IF (global .AND. (images_local .GT. 1)) THEN
     ! communicate all positions to all nodes (brute force, but simple)
#if defined(MPI) || defined (MPI_CHAIN)
        node = comm_chain%node_me

        R(:,:,:) = 0
        R(:,:,node) = posion
        CALLMPI_C( M_SUM_d( comm_chain, R(1,1,1), nions*3*images_local))

        F(:,:,:) = 0
        F(:,:,node) = force
        CALLMPI_C( M_SUM_d( comm_chain, F(1,1,1), nions*3*images_local))
#endif
      ELSE
        R(:,:,1) = posion
        F(:,:,1) = force
        IF(iu6>0) WRITE(iu6,*) "OPT: LBFGS force assigned"
      ENDIF
 
      IF (autoscale .EQV. .TRUE.) THEN

        ! Take a finite difference step down the force.
        IF (fdstep .EQV. .TRUE.) THEN
          direction = F/SQRT(SUM(F*F))

          IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, initial finite difference'

          Rold = R
          Fold = F

          ! Update the position.
          R = R + direction*finite_step

          fdstep = .FALSE.
        ELSE
          ! Compute invcurv.
          ! calculate curvature down direction
          direction = R - Rold
          DO im=1,images_local
            CALL set_pbc(direction(:,:,im))
          ENDDO
          direction_norm = SQRT(SUM(direction*direction))
          fp1 = SUM(Fold*direction/direction_norm)
          fp2 = SUM(F*direction/direction_norm)
          curvature = (fp1-fp2)/direction_norm

          ! "damping" makes invcurv smaller, which means that smaller
          ! steps will be taken along the descent direction. This conservative
          ! approach helps when the forces are noisy.
          invcurv = 1.0_q / curvature / damping
          IF(iu6>0) write(iu6,*) 'OPT: LBFGS, invcurv = ', invcurv

          IF (invcurv .LT. 0.0) THEN
            reset_flag = .TRUE.
            maxmove_flag = .TRUE.
            IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, reset due to negative curvature'
          ENDIF

          ! Reset the lbfgs memory
          IF(reset_flag) THEN
             itr = 0
             reset_flag = .FALSE.
          ELSE
             ! find new direction
             itr = itr + 1
             IF (itr .LE. memory) THEN
               IF(iu6>0) write(iu6,*) 'OPT: LBFGS, adding new vector'

               change_in_G(itr,:,:,:) = -(F-Fold)
               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, DeltaG:'
               IF (global .AND. (images_local .GT. 1)) THEN
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_G(itr,:,:,node)
               ELSE
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_G(itr,:,:,1)
               ENDIF

               change_in_R(itr,:,:,:) = R-Rold 
               ! apply periodic boundary conditions to these vector differences
               DO im=1,images_local
                 CALL set_pbc(change_in_R(itr,:,:,im))
               ENDDO

               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, DeltaR:'
               IF (global .AND. (images_local .GT. 1)) THEN
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_R(itr,:,:,node)
               ELSE
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_R(itr,:,:,1)
               ENDIF

               rho(itr) = 1.0_q/SUM(change_in_G(itr,:,:,:)*change_in_R(itr,:,:,:))
             ELSE
               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, shifting vectors'
               DO im = 1,memory-1
                 change_in_G(im,:,:,:) = change_in_G(im+1,:,:,:)
                 change_in_R(im,:,:,:) = change_in_R(im+1,:,:,:)
                 rho(im) = rho(im+1)
               ENDDO
               change_in_G(memory,:,:,:) = -(F-Fold)
               change_in_R(memory,:,:,:) = R-Rold
               ! apply periodic boundary conditions to these vector differences
               DO im=1,images_local
                 CALL set_pbc(change_in_R(memory,:,:,im))
               ENDDO
               rho(memory) = 1.0_q/SUM(change_in_G(memory,:,:,:)*change_in_R(memory,:,:,:))
             ENDIF
             IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, itr, rho: ',itr,rho(itr)
          ENDIF

          Rold = R
          Fold = F

          ! compute Ho*g
          IF (itr .LT. memory) THEN
            bound = itr
          ELSE
            bound = memory
          ENDIF
          IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, bound: ',bound
          direction = -1._q*F
          
          ! calculate Ho*g
          DO im=1,bound 
            jm = (bound + 1 - im)
            alpha(jm) = rho(jm)*SUM(change_in_R(jm,:,:,:)*direction)
            direction = direction - (alpha(jm)*change_in_G(jm,:,:,:))
          ENDDO
          direction = invcurv*direction !Ho=Ho*invcurv Ho is identity matrix
          DO im=1,bound
            beta = rho(im)*SUM(change_in_G(im,:,:,:)*direction)
            direction = direction + change_in_R(im,:,:,:)*(alpha(im)-beta)
          ENDDO

          ! direction down gradient
          direction = -1._q*direction

          step_size = ABS(SQRT(SUM(direction*direction)))
          IF (step_size .GT. SQRT(maxmove*maxmove*REAL(images_local))) THEN
            reset_flag = .TRUE.
            step_size = SQRT(REAL(images_local))*maxmove 
            ! If we make max move step take steepest descent step.
            direction = step_size * F/SQRT(SUM(F*F))
            IF(iu6>0) write(iu6,*) 'OPT: LBFGS, reset due to max move taking SD step'
          ENDIF

          IF (maxmove_flag .EQV. .TRUE.) THEN
            direction = F/SQRT(SUM(F*F))
            R = SQRT(REAL(images_local)) * maxmove * direction + R
            maxmove_flag = .FALSE.
          ELSE
            R = R + direction
          ENDIF

        ENDIF

      ! =========================
      ! MAIN ELSE DANO CODE BELOW
      ! =========================
      ELSE
        IF (fdstep) THEN  ! Finite difference step
          IF(iu6>0) WRITE(iu6,*) "OPT: LBFGS in fdstep"

          fdstep = .FALSE.
          optflag = .TRUE. ! keep control to do finite difference step

          ! check for reset of direction
          a1 = ABS(SUM(F*Fold))
          a2 = SUM(Fold*Fold)
  !        IF ((a1 .LE. 0.5_q*a2) .AND. (a2 .NE. 0.0_q)) THEN
          IF ((a1 .GT. 0.5_q*a2) .OR. (a2 .EQ. 0.0_q)) THEN
             reset_flag=.TRUE.
          ENDIF

          IF (lineopt .EQV. .FALSE.) reset_flag=.FALSE.
          IF (a2 .EQ. 0.0_q) reset_flag=.TRUE.

          ! Reset the lbfgs memory
          IF(reset_flag) THEN
             IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, init'
             itr=0
             reset_flag = .FALSE.
          ELSE
             ! find new direction
             itr = itr + 1
             IF (itr .LE. memory) THEN
               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, adding new vector'

               change_in_G(itr,:,:,:) = -(F-Fold)
               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, DeltaG:'
               IF (global .AND. (images_local .GT. 1)) THEN
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_G(itr,:,:,node)
               ELSE
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_G(itr,:,:,1)
               ENDIF

               change_in_R(itr,:,:,:) = R - Rold 
               ! apply periodic boundary conditions to these vector differences
               DO im=1,images_local
                 CALL set_pbc(change_in_R(itr,:,:,im))
               ENDDO

               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, DeltaR:'
               IF (global .AND. (images_local .GT. 1)) THEN
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_R(itr,:,:,node)
               ELSE
                 IF(iu6>0) WRITE(iu6,'(3f16.8)') change_in_R(itr,:,:,1)
               ENDIF

               rho(itr)=1.0_q/SUM(change_in_G(itr,:,:,:)*change_in_R(itr,:,:,:))
             ELSE
               IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, shifting vectors'
               DO im = 1,memory-1
                 change_in_G(im,:,:,:) = change_in_G(im+1,:,:,:)
                 change_in_R(im,:,:,:) = change_in_R(im+1,:,:,:)
                 rho(im) = rho(im+1)
               ENDDO
               change_in_G(memory,:,:,:) = -(F - Fold)
               change_in_R(memory,:,:,:) = R - Rold
               ! apply periodic boundary conditions to these vector differences
               DO im=1,images_local
                 CALL set_pbc(change_in_R(memory,:,:,im))
               ENDDO
               rho(memory) = 1.0_q/SUM(change_in_G(memory,:,:,:)*change_in_R(memory,:,:,:))
             ENDIF
             IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, itr, rho: ', itr, rho(itr)
          ENDIF

          Rold = R
          Fold = F

          ! compute Ho*g
          IF (itr .LT. memory) THEN
            bound = itr
          ELSE
            bound = memory
          ENDIF
          IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, bound: ',bound
          direction = -1._q*F
          
          ! calculate Ho*g
          IF (itr .NE. 0) THEN
            DO im=1,bound 
              jm = (bound + 1 - im)
              alpha(jm) = rho(jm)*SUM(change_in_R(jm,:,:,:)*direction)
              direction = direction - (alpha(jm)*change_in_G(jm,:,:,:))
            ENDDO
            direction = invcurv*direction !Ho=Ho*invcurv Ho is identity matrix
            DO im=1,bound
              beta = rho(im)*SUM(change_in_G(im,:,:,:)*direction)
              direction = direction + change_in_R(im,:,:,:)*(alpha(im)-beta)
            ENDDO
          ENDIF

          ! direction down gradient
          direction= -1._q*direction

          IF (lineopt) THEN ! line minimizer (using fdstep)
            direction = direction/SQRT(SUM(direction*direction))

            IF(iu6>0) WRITE(iu6,*) 'OPT: LBFGS, direction:'

            IF (global .AND. (images_local .GT. 1)) THEN
              IF(iu6>0) WRITE(iu6,'(3f16.8)') direction(:,:,node)
            ENDIF
            ! finite step down force
            R = R + direction*finite_step

          ELSE ! use hessian to make step
            step_size = ABS(SQRT(SUM(direction*direction)))
            IF (step_size .GT. SQRT(maxmove*maxmove*REAL(images_local))) THEN
              step_size = SQRT(maxmove*maxmove*REAL(images_local))
              direction = step_size*direction/SQRT(SUM(direction*direction))
            ENDIF
            R = R + direction
            fdstep = .TRUE.
            optflag = .FALSE. ! gives control back to method
          ENDIF


          IF(iu6>0) WRITE(iu6,*) "OPT: LBFGS, fdstep end"

        ELSE  ! Translation step for line-minimizer

          fdstep = .TRUE.
          optflag = .FALSE. ! gives control back to method

          ! calculate curvature down direction
          fp1 = SUM(Fold*direction)
          fp2 = SUM(F*direction)
          curvature = (fp1-fp2)/finite_step
          IF (curvature .LT. 0.0_q) THEN
            step_size = maxmove
          ELSE
            Favg = 0.5_q*(fp1 + fp2)
            step_size = Favg/curvature
            IF (ABS(step_size) .GT. maxmove) THEN
              step_size = SIGN(maxmove,step_size) - SIGN(finite_step,step_size)
            ELSE
              step_size = step_size - 0.5_q*finite_step  !(*)
            ENDIF
          ENDIF
          ! Move now from the configuration after the fd_step, so (*) has a "-" sign
          R = R + (direction*step_size)
        ENDIF
      ENDIF
      ! ========
      ! END DANO
      ! ========

      ! update posion
      IF (global .AND. (images_local .GT. 1)) THEN
        posion = R(:,:,node)
      ELSE
        posion = R(:,:,1)
      ENDIF

      IF(iu6>0) THEN
        WRITE(iu6,*)'OPT: LBFGS: positions in Cartesian coordinates' 
        WRITE(iu6,'(3f16.8)') posion
      ENDIF
      ! convert position back to direct cord
      CALL kardir(nions, posion, latt_b)
!#endif

      IF(iu6>0) WRITE(iu6,*) "OPT: LBFGS END subroutine"
    END SUBROUTINE lbfgs_step

!**********************************************************************
! lbfgs init
!**********************************************************************

    SUBROUTINE lbfgs_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,iu0
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions = T_INFO%nions
      iu0 = IO%IU0
      iu6 = IO%IU6

      ! read variables used for the LBFGS optimizer

      memory = 20
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'ILBFGSMEM','=','#',';','I', &
     &            memory,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''IBFGSMEM'' from file INCAR.'
         STOP
      ENDIF

      finite_step = 0.005_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'FDSTEP','=','#',';','F', &
     &            IDUM,finite_step,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''FDSTEP'' from file INCAR.'
         STOP
      ENDIF

      maxmove = 0.2_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'MAXMOVE','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''MAXMOVE'' from file INCAR.'
         STOP
      ENDIF
 
      global = .TRUE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LGLOBAL','=','#',';','L', &
     &            IDUM,RDUM,CDUM,global,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''LGLOBAL'' from file INCAR.'
         STOP
      ENDIF

      lineopt = .FALSE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LLINEOPT','=','#',';','L', &
     &            IDUM,RDUM,CDUM,lineopt,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''LLINEOPT'' from file INCAR.'
         STOP
      ENDIF

      autoscale = .TRUE.
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'LAUTOSCALE','=','#',';','L', &
     &            IDUM,RDUM,CDUM,autoscale,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''LAUTOSCALE'' from file INCAR.'
         STOP
      ENDIF

      invcurv = 0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'INVCURV','=','#',';','F', &
     &            IDUM,invcurv,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''INVCURV'' from file INCAR.'
         STOP
      ENDIF

      damping = 2.00_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'DAMPING','=','#',';','F', &
     &            IDUM,damping,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''DAMPING'' from file INCAR.'
         STOP
      ENDIF
 
      images_local = 0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'IMAGES','=','#',';','I', &
     &            images_local,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (iu0>=0) WRITE(iu0,*) 'Error reading item ''IMAGES'' from file INCAR.'
         STOP
      ENDIF

      IF (images_local .EQ. 0) images_local=1

      ! initialize the variables: vectors and matricies
      ALLOCATE(rho(memory),alpha(memory))
      IF (global .AND. (images_local .GT. 1)) THEN
        ALLOCATE(change_in_G(memory,3,nions,images_local))
        ALLOCATE(change_in_R(memory,3,nions,images_local))
        ALLOCATE(R(3,nions,images_local),Rold(3,nions,images_local))
        ALLOCATE(F(3,nions,images_local),Fold(3,nions,images_local))
        ALLOCATE(direction(3,nions,images_local))
      ELSE
        ALLOCATE(change_in_G(memory,3,nions,1))
        ALLOCATE(change_in_R(memory,3,nions,1))
        ALLOCATE(R(3,nions,1),Rold(3,nions,1))
        ALLOCATE(F(3,nions,1),Fold(3,nions,1),direction(3,nions,1))
      ENDIF
      change_in_G = 0._q
      change_in_R = 0._q
      Fold(:,:,1) = 0._q
      rho = 0._q
      alpha = 0._q
      itr = 0

      IF (iu6>=0) THEN
        WRITE(iu6,'(A5,A,I14.1)') 'OPT:',' LBFGS, MEMORY ',memory
        WRITE(iu6,'(A5,A,I14.1)') 'OPT:',' LBFGS, IMAGES ',images_local
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' LBFGS, MAXMOVE',maxmove
        WRITE(iu6,'(A5,A,L7)') 'OPT:',' LBFGS, LGLOBAL',global
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' LBFGS, INVCURV',invcurv
        WRITE(iu6,'(A5,A,L7)') 'OPT:',' LBFGS, LLINEOPT',lineopt
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' LBFGS, FDSTEP ',finite_step
        WRITE(iu6,'(A5,A,L7)') 'OPT:',' LBFGS, LAUTOSCALE',autoscale
        WRITE(iu6,'(A5,A,F14.6)') 'OPT:',' LBFGS, DAMPING',damping
      END IF

      reset_flag = .TRUE.
      fdstep = .TRUE.

    END SUBROUTINE lbfgs_init

!**********************************************************************
! Vector Functions
!**********************************************************************

!======================================================================
! Sets a vector to have the smallest length consistent the the periodic
! boundary conditions.
!======================================================================
      SUBROUTINE set_pbc(v)
      REAL(q) :: v(3,nions)
      CALL kardir(nions,v,car2dir)
      v = MOD(v + 100.5_q, 1._q) - 0.5_q
      CALL dirkar(nions,v,dir2car)
      END SUBROUTINE set_pbc

  END MODULE lbfgs

