#include "symbol.inc"
!**********************************************************************
!
! Module for running dynamics on system
!
! Version 1.01, July 2006
!
! TO DO
!
!**********************************************************************

  MODULE dynamic
    USE prec
    USE lattice

    IMPLICIT NONE
    PRIVATE :: bivarate
    PUBLIC :: dynamic_step, dynamic_init  !CALL opt_init from chain_init

    INTEGER :: nions,iu6,thermosdat
    REAL(q),Alocateable  R(:,:)
    REAL(q) :: dt,temp

    ! varables for Langevin dynamics 
    REAL(q) :: c0,c1,c2
    REAL(q),ALLOCATABLE :: Fold(:,:),velocity(:,:),randR(:,:),randV(:,:)


!**********************************************************************
!
! dynamics rutine using verlet velocity style dynamics
!
!**********************************************************************
  CONTAINS
    SUBROUTINE dynamic_step(optflag,tinfo_nions,posion,toten,force, &
                        latt_a,latt_b,tinfo_iu6)

      INTEGER :: tinfo_nions,tinfo_iu6
      REAL(q) :: posion(3,nions),toten,force(3,nions)
      REAL(q) :: latt_a(3,3),latt_b(3,3)
      LOGICAL optflag

      REAL(q) :: a1,a2,curvature

      R = posion

      ! Convert the position into Cartesian coordinates
      CALL dirkar(nions,R,latt_a)

      ! give contral back to method 
      optflag = .FALSE.

      IF (thermosdat .EQ. 1) THEN
        posion = posion + c1*dt*velocity + c2*dt*dt*Fold/mass + randR
        velocity = c0*velocity+(c1-c2)*dt*Fold/mass+c2*dt*force/mass+randV
        Fold = force
      ELSE


      ! convert position back to direct cord
      CALL kardir(nions,R,latt_b)
      ! update posion
      posion = R

    END SUBROUTINE dynmaic_step

!**********************************************************************
! bivariate
!
!  returns 2 random numbers from bivariate Gaussian dist. with specified
!  std's and correlation coef.
!  i.e., 
! p(g1,g2)/norm 
! = exp(-(sig2**2*g1**2+sig1**2*g2**2)/(2[sig1**2+sig2**2-(c12*sig1*sig2)**2])
! = exp(-(var2*g1**2+var1*g2**2)/(2(var1+var2-c12**2*var1*var2))
! e.g.,  see p 28 of Wax (Chandrasekhar article, Eq.(178) )
! uses Box-Muller transformation;  see p. 203 of Num. Rec., 1st ed.
! and then covariance formula from App. G of Allen and Tildesley, p. 348.
!**********************************************************************
    SUBROUTINE bivariate(sig1,sig2,c12,g1,g2)

      REAL(q) :: sig1,sig2,c12,g1,g2

1      v1=2.0_q*rand()-1.0_q   ! pick coords in square (-1,1)x(-1,1)
      v2=2.0_q*rand()-1.0_q
      v=v1*v1+v2*v2
      IF(v.GE.1.0_q) GOTO 1     ! discard pairs that are not in unit circle
      fac=SQRT(-2.0_q*LOG(v)/v)  ! Box-Muller transformation
      gset1=v1*fac      
      gset2=v2*fac     

      g1=sig1*gset1
      g2=sig2*(c12*gset1 + SQRT(1.0_q-c12*c12)*gset2)
      RETURN
     END SUBROUTINE  


!**********************************************************************
!
! dynamic initilizer
!
!**********************************************************************
    SUBROUTINE dynamic_init(T_INFO,IO)
      USE base
      USE poscar
      TYPE(in_struct) :: IO
      TYPE(type_info) :: T_INFO


      INTEGER IDUM,IERR,N
      CHARACTER*1 CHARAC
      COMPLEX(q) CDUM
      LOGICAL LDUM
      REAL(q) RDUM

      nions=T_INFO%nions
      iu6=IO%IU6

      REAL(q) :: gdt

      ! read in varables used for dynamics

      thermosdat=0
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'THERMOSDAT','=','#',';','I', &
     &            thermosdat,RDUM,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''THERMOSDAT'' from file INCAR.'
         STOP
      ENDIF

      dt = 0.01_q
      CALL RDATAB(.TRUE.,'INCAR',IO%IU5,'TIMESTEP','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''TIMESTEP'' from file INCAR.'
         STOP
      ENDIF

      temp = 300.0_q
      CALL RDATAB(.TRUE.,0'INCAR',IO%IU5,'MDTEMP','=','#',';','F', &
     &            IDUM,temp,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''MDTEMP'' from file INCAR.'
         STOP
      ENDIF

      gam = 0.1_q
      CALL RDATAB(.TRUE.,0'INCAR',IO%IU5,'MDGAMMA','=','#',';','F', &
     &            IDUm,gam,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''MDGAMMA'' from file INCAR.'
         STOP
      ENDIF

      alpha = 0.5_q
      CALL RDATAB(.TRUE.,0'INCAR',IO%IU5,'MDALPHA','=','#',';','F', &
     &            IDUM,alpha,CDUM,LDUM,CHARAC,N,1,IERR)
      IF (((IERR/=0).AND.(IERR/=3)).OR. ((IERR==0).AND.(N<1))) THEN
         IF (IO%IU0>=0) &
              WRITE(iu6,*)'Error reading item ''MDALPHA'' from file INCAR.'
         STOP
      ENDIF

      ! set up thermosdat parameters
      IF (thermosdat .EQ. 0) .OR. (gam .EQ. 0.0_q) THEN
        c0 = 1.0_q
        c1 = 1.0_q
        c2 = 0.5_q
      ELSE 
        !coeficients for Brownian integration (Allen and Tildesley, p. 261)
        gdt = mdgamma*dt
        c0 = EXP(-gdt)
        c1 = (1.0_q-c0)/(gdt)
        c2 = (1.0_q-c1)/(gdt)
        sig2r = temp*k/mass/mdgamma/mdgamma *
     x      (2.0_q*gdt-3.0_q+4.0_q*EXP(-gdt)-EXP(-2.0_q*gdt))
        sig2v=temp*k/mass*(1.0_q-EXP(-2.0_q*gdt))
        sigrsigvcrv=temp*k/mass/mdgamma*(1.0_q-EXP(-gdt))**2
        sigr=SQRT(sig2r)
        sigv=SQRT(sig2v)
        crv=sigrsigvcrv/(sigr*sigv)
      ENDIF
      !variances and correlation coef. for x and v displacements (A&T, p. 262)

      ! initialize the variables: vectors and matricies
      ALLOCATE(Fold(3,nions),velocity(3,nions),direction(3,nions))


    END SUBROUTINE dynamic_init
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
      v1=MOD(v1+100.5_q,1._q)-0.5_q
      CALL dirkar(nions,v1,dir2car)
      END SUBROUTINE set_pbc
!======================================================================
! Returns a unit vector along v1
!======================================================================
      FUNCTION return_unit(V1)
      real(q) :: v1(3,nions)
      real(q),dimension(3,nions) :: return_unit
      return_unit=v1*(1._q/SQRT(SUM(v1*v1)))
      END FUNCTION return_unit
!======================================================================
! Sets V1 to be a unit vector
!======================================================================
      SUBROUTINE set_unit(V1)
      REAL(q) :: v1(3,nions)
      v1=return_unit(v1)
      END SUBROUTINE set_unit
!======================================================================
! Vector projection of v1 on v2
!======================================================================
      FUNCTION vproj(v1,v2)
      REAL(q) :: v1(3,nions),v2(3,nions),vproj(3,nions)
      vproj=v2*SUM(v1*v2)/SQRT(SUM(v2*v2))
      END FUNCTION vproj


  END MODULE dynamic


