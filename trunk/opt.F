#include "symbol.inc"
!**********************************************************************
! RCS:  $Id: opt.F,v 1.5 2006-06-26 23:51:08 graeme Exp $
!
! Module which implements our optimizers
!
! Version 1.01, June 2006
!
!**********************************************************************

  MODULE opt
    USE prec
    USE lattice

    IMPLICIT NONE
    private :: sd_step, cg_step, qm_step, lbfgs_step, bfgs_step
    private :: sd_init, cg_init, qm_init, lbfgs_init, bfgs_init

    public :: opt_step, opt_init  !call opt_init from chain_init

    INTEGER :: nions,iu6
    REAL(q),ALLOCATABLE :: R(:,:),F(:,:)
    REAL(q),ALLOCATABLE :: Rold(:,:),Fold(:,:)

    REAL(q) :: dir2car(3,3),car2dir(3,3)
    REAL(q) :: finite_step ! ...
    LOGICAL :: fdstep, optflag_local

!**********************************************************************
!
! Routine to optimize
!
!**********************************************************************

  CONTAINS
    SUBROUTINE opt_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)
      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      ! local variables

      nions=tinfo_nions
      iu6=tinfo_iu6
      dir2car=latt_a
      car2dir=latt_b
      optflag_local=optflag

      ! if optflag is false, do nothing
      if(.not. optflag) return

      ! things we want to do for all optimizers?

     if (iopt==1) then
        call lbfgs_step()
      elseif (iopt==2) then
        call cg_step()
      elseif (iopt==3) then
        call qm_step()
      elseif (iopt==4) then
        call sd_step()
      elseif (iopt==5) then
        call bfgs_step()
      endif

      ! call periodic boundary conditions
      optflag=optflag_local

    END SUBROUTINE opt_step


    SUBROUTINE opt_init(T_INFO,IO)
      type(in_struct) :: IO
      type(type_info) :: T_INFO

      INTEGER IDUM,IERR,N,ibrion
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM,potim

      nions=T_INFO%nions
      iu6=IO%IU6

      ! read iopt (0=use vasp optimizers, 1=lbfgs, 2=cg, 3=qm, 4=bfgs, 5=sd)

      iopt=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'IOPT','=','#',';','F', &
     &            IOPT,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''IOPT'' from file INCAR.'
         STOP
      ENDIF
      ! write out the optimizer chosen

      if (iopt==0) return

      potim=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'potim','=','#',';','F', &
     &            IDUM,potim,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''POTIM'' from file INCAR.'
         STOP
      ENDIF

      ibrion=0
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'IBRION','=','#',';','F', &
     &            IBRION,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''IBRION'' from file INCAR.'
         STOP
      ENDIF

      IF(potim/=0 .OR. ibrion/=3) THEN
        WRITE(IO%IU6,*) 'OPT: Must set IBRION=3 and POTIM=0 for IOPT>0'
        STOP
      ENDIF

      ! read optimization parameters

      finite_step=0.001_q
      call RDATAB(.True.,'INCAR',IO%IU5,'RSTEP','=','#',';','I', &
     &            IDUM,finite_step,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR.((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(IO%IU0,*)'Error reading item ''RSTEP'' from file INCAR.'
         STOP
      ENDIF

      ! allocate memory for local F and R

      allocate(F(3,Nions),R(3,Nions))
      allocate(Fold(3,Nions),Rold(3,Nions))

      ! optimizer specific initialization
      
      if (iopt==1) then
        call lbfgs_init()
      elseif (iopt==2) then
        call cg_init()
      elseif (iopt==3) then
        call qm_init()
      elseif (iopt==4) then
        call sd_init()
      elseif (iopt==5) then
        call bfgs_init()
      endif

    end subroutine opt_init

!**********************************************************************
!
! Limited-memory bfgs method
!
!**********************************************************************

    subroutine lbfgs_init()

      REAL(q),ALLOCATABLE,DIMENSION(:,:,:),SAVE :: change_in_force,change_in_posion

      INTEGER IDUM,IERR,N
      INTEGER,SAVE memory
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      memory=4
      CALL RDATAB(IO%LOPEN,'INCAR',IO%IU5,'NLBFGSMEM','=','#',';','F', &
     &            memory,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(TIU0,*)'Error reading item ''NBFGSMEM'' from file INCAR.'
         STOP
      ENDIF

      ! initialize the variables: vectors and matricies
      allocate(change_in_force(memory,3,Nions),change_in_posion(memory,3,Nions))

      fdstep=.TRUE.

    end subroutine lbfgs_init

    subroutine lbfgs_step()
      REAL(q) :: change_in_force,change_in_posion !matricies

      REAL(q) :: rho,q

      REAL(q) :: alpha,beta ! these are vectors
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,avg_force
      REAL(q) :: finite_step,step_size
      INTEGER :: itr,incr,bound,im,jm

      if (finite_difference_step_flag) then
        finite_difference_step_flag=.false.
        ! check for reset of direction
        a1 = abs(sum(force*force_old))
        a2 = sum(force_old*force_old)
        if (a1.lt.0.5_q*a2) then
           ! find new direct
           change_in_force(itr,:,:) = F-Fold
           change_in_posion(itr,:,:) = R-Rold  ! need boundary conditions
           rho(itr) = 1.0_q/sum(change_in_force(itr,:,:)*change_in_posion(itr,:,:))
           itr = itr+1
        else 
           ! reset down the force
           itr = 1 
           ! ? reset the change_in_pos/force vectors?
        endif 
        
        posion_old = posion
        force_old = force
        ! compute Ho*g
        if (itr.lt.memory) then 
          incr = 0
          bound = itr
        else 
          incr = itr-memory
          bound = memory
        endif
        q = force
        do im=0,bound-1  ! change the range to 1..N
          jm = (bound-im)+incr
          alpha(bound-im) = rho(jm)*sum(change_in_posion(jm)*q)
          q = q-alpha(bound-im)*change_in_force(jm)
        enddo
        direction = Ho*q
        do im=0,bound-1
          jm = im+incr
          beta = rho(jm)*sum(change_in_force(jm)*direction
          direction = direction + change_in_posion*(alpha(im)-beta)
        enddo 
        direction = return_unit(direction)

        ! finite step down force
        posion = posion+direction*step_size
      else  
        finite_difference_step_flag=.true.
        ! calculate curvature down direction
        fp1 = sum(force_old*direction)
        fp2 = sum(force*direction)
        avg_force = (fp1+fp2)/2.0_q
        curvature = (fp1-fp2)/finite_step
        if (curvature<0.0_q) then 
          step_size=maxmove
        else
          step_size = avg_force/curvature
          if(abs(step_size)>maxmove) then
            step_size=sign(maxmove,step_size)
          else
            step_size = step_size-(finite_step/2.0)
          endif
        endif
        posion = posion+(direction*step_size)
      endif

    return

    end subroutine lbfgs_step


!**********************************************************************
!
! Conjugate gradient method
!
!**********************************************************************

    subroutine cg_init()
    ! load the variables: dR, maxmove
    ! initialize flag: finite_diff_step_flag
    end subroutine cg_init

    subroutine cg_step(posion,toten,force)
!      LOGICAL,INTENT(INOUT) :: fd_step
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q),PARAMETER :: gamma=0.5_q
      REAL(q) :: c1,c2,f1d,f2d,fr,cr,step,s

!      IF (fd_step) THEN
      IF (finite_difference_step_flag) THEN
        finite_difference_step_flag=.false.
        f1=Feff
        c1=ABS(SUM(f1*fco))
        c2=SUM(fco*fco)
        IF (c1 < gamma*c2) THEN
          s=SUM(f1*(f1-fco))/c2
        ELSE
          s=0.0_q
        END IF
        dc=f1+s*dco
        dcu=dc/SQRT(SUM(dc*dc))
        fco=f1
        dco=dc
! Move from the original configuration
        R=R0+dcu*dR
      ELSE
!        fd_step=.true.
        finite_difference_step_flag=.true.
        f2=Feff
        f1d=SUM(F1*dcu)
        f2d=SUM(f2*dcu)
        cr=(f1d-f2d)/dR
        IF (cr < 0.0_q) THEN
         step=maxmove
        ELSE
          fr=0.5_q*(f1d+f2d)
          step=fr/cr
          IF (ABS(step) > maxmove) THEN
            step=SIGN(maxmove,step)
          ELSE
            step=step-0.5_q*dR     ! (*)
          END IF
        END IF
! Move now from the configuration after the fd_step, so (*) has a "-" sign
        R=R+dcu*step
      END IF

    RETURN


    end subroutine cg_step


!**********************************************************************
!
! Qucik-min method
!
!**********************************************************************

    subroutine qm_init()
    ! load the variables: dt, maxmove
    ! set velocity
    end subroutine qm_init

    subroutine qm_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: velocity,step

      if ((sum(force*velocity)).ge.0.0_q) then
        velocity = vproj(velocity,force)
      else
        velocity = velocity*0.0_q
      endif

      ! Euler step
      velocity = velocity*dt*force/mass
      step = dt*velocity
      if (sqrt(sum(step*step)).ge.maxmove) then
        step = maxmove*return_unit(step)
      endif

      posion = posion+step

    return ! do i need this?
    end subroutine qm_step


!**********************************************************************
!
! Steepest-Decent method
!
!**********************************************************************

    subroutine sd_init()
    ! load the variables: alpha, maxmove
    end subroutine sd_init

    subroutine sd_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: step
    
      step = force * alpha
      if (sqrt(sum(step*step)).ge.maxmove then
        step = maxmove*return_unit(step)
      endif

      posion = posion+step

    return ! do i need this?
    end subroutine sd_step

!**********************************************************************
!
! full bfgs method
!
!**********************************************************************

    subroutine bfgs_init()
    ! load the variables: dR, maxmove, memory
    ! initialize flag: finite_diff_step_flag
    end subroutine bfgs_init

    subroutine bfgs_step(posion,toten,force)
      REAL(q) :: posion(3,nions),toten,force(3,nions)

      REAL(q) :: posion_old,force_old
      REAL(q) :: a1,a2,direction,curvature
      REAL(q) :: fp1,fp2,avg_force
      REAL(q) :: finite_step,step_size

      if (finite_difference_step_flag) then
        finite_difference_step_flag=.false.
        ! check for reset of direction
        a1 = abs(sum(force,force_old))
        a2 = sum(force_old,force_old)
        if (a1.lt.0.5_q*a2) then
           ! find new direct
        else 
           ! the hessian to Idenity
        endif 
        
        posion_old = posion
        force_old = force
        ! compute Ho*g


        direction = return_unit(direction)

        ! finite step down force
        posion = posion+direction*step_size
      else  
        finite_difference_step_flag=.true.
        ! calculate curvature down direction
        fp1 = sum(force_old*direction)
        fp2 = sum(force*direction)
        avg_force = (fp1+fp2)/2.0_q
        curvature = (fp1-fp2)/finite_step
        if (curvature<0.0_q) then 
          step_size=maxmove
        else
          step_size = avg_force/curvature
          if(abs(step_size)>maxmove) then
            step_size=sign(maxmove,step_size)
          else
            step_size = step_size-(finite_step/2.0)
          endif
        endif
        posion = posion+(direction*step_size)
        ! update the Hessian
 
      endif

    return


    end subroutine bfgs_step



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

END MODULE opt
