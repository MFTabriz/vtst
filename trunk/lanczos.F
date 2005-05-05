#include "symbol.inc"
!**********************************************************************
! RCS:   $Id: lanczos.F , v0.01 June 3th 2004
!
! This module implements the Lanczos method for saddle point finding.
!
! More work is needed, there WILL be A LOT of bugs. Output needs to be
! added and refined. 
! The module was developed on an IBM AIX system so expect problems when
! porting to other systems. Both serial and parallel versions should
! work on AIX machines.
!
! Andri Arnaldsson
! andri@u.washington.edu
!
!**********************************************************************!

  MODULE lanczos
    USE prec
    USE main_mpi
    USE poscar
    USE lattice

    IMPLICIT NONE
    SAVE
    PRIVATE
    PUBLIC :: lanczos_force,lanczos_init

    TYPE(in_struct) :: lio
    TYPE(type_info) :: linfo
    INTEGER :: nions,iu6,iu0,nl
    INTEGER,PARAMETER :: lanout=831
    REAL(q),ALLOCATABLE,DIMENSION(:) :: d,e,aa,bb
    REAL(q),ALLOCATABLE,DIMENSION(:,:) :: F,R,w,qq,qqold,F0,R0,z,Fpar,Feff,vel,dc,dco, &
  &                                       fco,f1,f2,dcu
    REAL(q),ALLOCATABLE,DIMENSION(:,:,:) :: PP
    REAL(q) :: dR,ltol,dt,maxmove
    LOGICAL :: ifcg

    CONTAINS
!--------------------------------------------------------------------------------------!

    SUBROUTINE lanczos_force(nions,posion,toten,tifor,latt_a,latt_b,iu6)
      INTEGER :: nions,iu6
      REAL(q),DIMENSION(3,nions) :: posion,tifor 
      REAL(q),DIMENSION(3,3) :: latt_a,latt_b
      REAL(q) :: toten
      REAL(q),EXTERNAL :: rane

      REAL(q),DIMENSION(3,3) :: A,B
      REAL(q),SAVE :: alpha,beta,eigold,eig
      REAL(q) :: U
      INTEGER,SAVE :: it=0,itr=0
      INTEGER :: i,j
      LOGICAL,SAVE :: first=.true.,new,iterate,converged,fd_step
      LOGICAL :: ertil

      A=latt_a
      B=latt_b
      U=toten
      F=tifor
      R=posion
      CALL dirkar(nions,R,A)

!!!!!!!!!!!! TEST CG !!!!!!!!!!!!!!!!!!
!      if (first) then                !
!        first=.false.                !
!        fd_step=.true.               !
!      end if                         !
!      Feff=F                         !
!      R0=R                           !
!      call cg(fd_step)               !
!      posion=R                       !
!      call kardir(nions,posion,B)    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!! TEST QUICK-MIN !!!!!!!!!!!
!      Feff=F                         !
!      R0=R                           !
!      call quickmin()                !
!      posion=R                       !
!      call kardir(nions,posion,B)    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (first) THEN
        first=.false.
        new=.true.
! Open the ouput file for the run
!        OPEN(lanout,FILE='lanczos.out',ACTION='write',STATUS='replace')        

        OPEN(lanout,FILE=DIR_APP(1:DIR_LEN)//'lanczos.out',STATUS='unknown')

        IF (iu6 > 0) CALL output(0)
        IF (iu6 > 0) CALL output(50)
        IF (iu6 > 0) CALL output(80)
! If the <MODECAR> exists then read in the initial guess for the mode
! otherwise make a random guess.
        IF (iu6 > 0) INQUIRE(FILE='MODECAR',EXIST=ertil)
        IF (ertil) THEN
          IF (iu6 > 0) THEN
            OPEN(210,FILE='MODECAR',ACTION='read',STATUS='old')
            READ(210,'(3ES20.10)') (w(:,i) , i=1,nions)
            CLOSE(210)
          END IF
        ELSE
          DO i=1,3
            DO j=1,nions
              IF (linfo%lsfor(i,j)) THEN
                w(i,j)=rane()-0.5_q
              END IF
            END DO
          END DO
        END IF
      END IF

      IF (new) THEN
        new=.false.
        iterate=.true.
        converged=.false.
        IF (ifcg) fd_step=.true.
        beta=SQRT(SUM(w*w))        
! Save the forces and coordinates at the point we are in.
        F0=F
        R0=R
        itr=itr+1
        IF (iu6 > 0) CALL output(100,itr)
! Leave to calculate new force values
        qqold=0.0_q
        qq=w/beta
        it=1  
        PP(:,:,it)=qq
        R=R0+qq*dR   
      ELSE
        IF (iterate) THEN
          z=F-F0
          w=z-beta*qqold
          alpha=SUM(w*qq)
          d(it)=alpha
          w=w-alpha*qq
          beta=SQRT(SUM(w*w))
          e(it)=beta 
! Check the eigenvalues          
          IF (it > 1) THEN
            aa(1:it)=-d(1:it)/dR
            bb(1:it)=-e(1:it)/dR
            CALL eigenvalue(aa(1:it),bb(1:it),it)
            CALL sort(it,aa(1:it),bb(1:it))
            eig=aa(1)
            converged=(ABS((eig-eigold)/eigold) < ltol)
            IF (iu6 > 0) CALL output(200,itr,it,eig,eigold)
          END IF             
        END IF
! If the size of the Lanczos matrix has reached its maximum allowed dimension (nl),
! the continue on with the latest (unconverged) eigenvalue and print a warning 
! message to "lanczos.out".
        IF ((it == nl) .AND. (.NOT. converged)) THEN
          converged=.true.
          IF (iu6 > 0) CALL output(250,itr)
        END IF

        IF (converged) THEN
! Only need to find the mode ones but goes through here twice if CG is used, so ...
          IF (iterate) THEN
            CALL eigenvector(it,eig)
            IF (iu6 > 0) CALL output(300,itr)
          END IF        
          IF (ifcg) THEN
            iterate=.false.
            IF (.NOT. fd_step) F0=F   ! for feffective 
            CALL feffective(eig)
            CALL cg(fd_step)
            IF ((.NOT. fd_step) .AND. (iu6 > 0)) CALL output(400,itr,it,eig,U=U)   ! only output from the point, not after fd step
            IF (fd_step) new=.true.
          ELSE 
            CALL feffective(eig)
            IF (iu6 > 0) CALL output(400,itr,it,eig,U=U)
            CALL quickmin()
! If the eigenvalue has converged, exit to a new step
            new=.true.
          END IF
        ELSE
          IF (it > 1) THEN
            eigold=eig
          ELSE
            eigold=-alpha/dR
          END IF
          it=it+1
! Leave to calculate new force values
          qqold=qq
          qq=w/beta
          PP(:,:,it)=qq
          R=R0+qq*dR
        END IF
      END IF                       ! end IF (new)

      posion=R
      CALL kardir(nions,posion,B)
  
      RETURN
    END SUBROUTINE lanczos_force

!--------------------------------------------------------------------------------------!

    SUBROUTINE lanczos_init(t_info,io)
      TYPE(in_struct) :: io
      TYPE(type_info) :: t_info

      lio=io
      linfo=t_info
      iu6=lio%iu6
      iu0=lio%iu0
      nions=t_info%nions

      CALL read_variables(lio)

      ALLOCATE(F(3,nions),R(3,nions),w(3,nions),qq(3,nions),qqold(3,nions),F0(3,nions),&
  &               R0(3,nions),z(3,nions),Feff(3,nions))
      ALLOCATE(PP(3,nions,nl))
      ALLOCATE(d(nl),e(nl),aa(nl),bb(nl))

      F=0.0_q
      R=0.0_q
      w=0.0_q
      qq=0.0_q
      qqold=0.0_q
      PP=0.0_q
      F0=0.0_q
      R0=0.0_q
      z=0.0_q
      d=0.0_q
      e=0.0_q
      Feff=0.0_q
      IF (ifcg) THEN
        ALLOCATE(dc(3,nions),dco(3,nions),fco(3,nions),f1(3,nions),f2(3,nions),        &    
  &              dcu(3,nions))
        dc=0.0_q
        dco=0.0_q
        fco=0.0_q
        f1=0.0_q 
        f2=0.0_q   
        dcu=0.0_q
      ELSE
        ALLOCATE(vel(3,nions))
        vel=0.0_q
      END IF

    RETURN
    END SUBROUTINE lanczos_init

!--------------------------------------------------------------------------------------!

    SUBROUTINE feffective(eig)
      REAL(q),INTENT(IN) :: eig
      
      REAL(q),DIMENSION(3,nions) :: Fpar 

      Fpar=w*SUM(w*F0)
      IF (eig < 0.0_q) THEN
        Feff=F0-2.0_q*Fpar
      ELSE
        Feff=-Fpar
      END IF

    RETURN
    END SUBROUTINE feffective

!--------------------------------------------------------------------------------------!

    SUBROUTINE cg(fd_step)
      LOGICAL,INTENT(INOUT) :: fd_step

      REAL(q),PARAMETER :: gamma=0.5_q
      REAL(q) :: c1,c2,f1d,f2d,fr,cr,step,s

      IF (fd_step) THEN
        fd_step=.false.
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
        fd_step=.true.
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
    END SUBROUTINE cg

!--------------------------------------------------------------------------------------!
        
    SUBROUTINE quickmin()
        
      REAL(q),DIMENSION(3,nions) :: dv,dRq
      REAL(q) :: t
        
      dv=Feff*dt
      WHERE(Feff*vel < 0.0_q) vel=0.0_q
      vel=Feff*SUM(vel*Feff)/SUM(Feff*Feff)
      vel=vel+dv
      dRq=vel*dt
      t=SQRT(SUM(dRq*dRq)) 
      IF (t > MaxMove) dRq=dRq/t*MaxMove
      R=R0+dRq
        
    RETURN
    END SUBROUTINE quickmin
        
!------------------------------------------------------------------------------------!
     
    SUBROUTINE read_variables(io)
      TYPE(in_struct) :: io
      integer :: IDUM,IERR,Nint
      character*1 :: CHARAC
      complex(q) :: CDUM
      logical :: LDUM
      real(q) :: RDUM
      
      ltol=1.0e-2_q
      call RDATAB(.True.,'INCAR',IO%IU5,'Sltol','=','#',';','F', &
     &            IDUM,ltol,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dR=1.0e-3_q
      call RDATAB(.True.,'INCAR',IO%IU5,'SdR','=','#',';','F', &
     &            IDUM,dR,CDUM,LDUM,CHARAC,Nint,1,IERR)
      maxmove=0.3_q
      call RDATAB(.True.,'INCAR',IO%IU5,'Smaxmove','=','#',';','F', &
     &            IDUM,maxmove,CDUM,LDUM,CHARAC,Nint,1,IERR)
      dt=0.15_q  
      call RDATAB(.True.,'INCAR',IO%IU5,'Sdt','=','#',';','F', &
     &            IDUM,dt,CDUM,LDUM,CHARAC,Nint,1,IERR)
      nl=20 
      call RDATAB(.True.,'INCAR',IO%IU5,'Snl','=','#',';','I', &
     &            nl,RDUM,CDUM,LDUM,CHARAC,Nint,1,IERR)
      ifcg=.true.   
      call RDATAB(.True.,'INCAR',IO%IU5,'ifcg','=','#',';','L', &
     &            IDUM,RDUM,CDUM,ifcg,CHARAC,Nint,1,IERR)


!      if (iu6 > 0) write(*,*) ltol
!      if (iu6 > 0) write(*,*) dR
!      if (iu6 > 0) write(*,*) maxmove
!      if (iu6 > 0) write(*,*) dt
!      if (iu6 > 0) write(*,*) nl
!      if (iu6 > 0) write(*,*) ifcg

    RETURN        
    END SUBROUTINE read_variables

!------------------------------------------------------------------------------------!

    SUBROUTINE output(line,itr,it,eig,eigold,U)
      INTEGER,INTENT(IN) :: line
      INTEGER,OPTIONAL,INTENT(IN) :: itr,it
      REAL(q),OPTIONAL,INTENT(IN) :: eig,eigold,U

      INTEGER :: i

      SELECT CASE (line)
        CASE (0)
          WRITE(lanout,'(A)') 'Echo control variables:'
          WRITE(lanout,'(A,1I3)') 'Maximum size of Lanczos matrix .............. nl=',nl
          WRITE(lanout,'(A,1ES8.1)') 'Finite difference step length ............... dR=',dR
          WRITE(lanout,'(A,1ES8.1)') 'Tolerance for eigenvalue convergence ...... ltol=',ltol
          WRITE(lanout,'(A,1ES8.1)') 'Time step length for Quick-Min .............. dt=',dt
          WRITE(lanout,'(A,1ES8.1)') 'Maximum total movement in one step ..... maxmove=',maxmove
          WRITE(lanout,'(A,1L2)') 'Use conjugate gradients when minimizing ... ifcg=',ifcg
          WRITE(lanout,*) ' '
        CASE (50)
          WRITE(lanout,'(A)') 'eig   Step#   Iteration     Eig           EigOld       |(Eig-EigOld)/EigOld|'
          WRITE(lanout,'(A)') 'eig  -------------------------------------------------------------------------'
        CASE (80)
          WRITE(lanout,'(A)') 'conv  Step#    Iteration      Energy         Max|Feff|      Eig'
          WRITE(lanout,'(A)') 'conv ----------------------------------------------------------'
        CASE (100)
          WRITE(lanout,'(/,A,I4,A)') 'Point ',itr,':'
          WRITE(lanout,'(A)') '------------------------------------------------------'
          WRITE(lanout,'(/,A)') ' Coordinates:'
          WRITE(lanout,'(1A5,1I4,1F16.12,2F20.12)') ('Coo  ',itr,R(:,i) , i=1,nions)
          WRITE(lanout,'(/,A)') ' Forces:'
          WRITE(lanout,'(1A5,1I4,3ES22.12)') ('For  ',itr,F(:,i) , i=1,nions)          
          WRITE(lanout,*) ' '
        CASE (200)
          WRITE(lanout,'(1A5,1I4,6X,1I4,3X,3G16.8)')  &
  &                   'eig  ',itr,it,eig,eigold,ABS((eig-eigold)/eigold)
        CASE (250)
          WRITE(lanout,'(A,A,1I4)') 'eig  ','Warning !!! Unconverged eigenvalue at point ',itr 
        CASE (300)
          WRITE(lanout,'(/,A)') ' Lowest Mode'
          WRITE(lanout,'(1A5,1I4,3ES20.10)') ('Low  ',itr,w(:,i) , i=1,nions)
          WRITE(lanout,*) ' '
        CASE (400)
          WRITE(lanout,'(1A6,1I4,6X,1I4,3X,1F14.5,4X,1ES12.2,2X,1F12.5)')      &
    &             'conv  ',itr,it,U,MAXVAL(ABS(Feff)),eig
      END SELECT
 
    RETURN
    END SUBROUTINE output

!------------------------------------------------------------------------------------!

    SUBROUTINE eigenvector(il,eig)
      INTEGER,INTENT(IN) :: il
      REAL(q),INTENT(IN) :: eig
      REAL(q),EXTERNAL :: rane
 
      INTEGER :: i
      REAL(q) :: shift,t
      LOGICAL :: first_ii

! Eigevalue is converged, find the eigenvector using an inverse iteration scheme
      d=-d/dR
      e=-e/dR 
      shift=-eig+0.0001_q  
      d=d+shift
! Cholesky factorize the matrix
      DO i=2,il
        t=e(i-1)
        e(i-1)=t/d(i-1)
        d(i)=d(i)-t*e(i-1)
      END DO
      first_ii=.true.
! Inverse interation
      DO 
        IF (first_ii) THEN 
          first_ii=.false.
          DO i=1,il
            bb(i)=rane()-0.5_q
          END DO
          bb(1:il)=bb(1:il)/SQRT(DOT_PRODUCT(bb(1:il),bb(1:il)))
          aa(1:il)=bb(1:il)
        END IF
        DO i=2,il
          bb(i)=bb(i)-e(i-1)*bb(i-1)
        END DO   
        bb(il)=bb(il)/d(il)
        DO i=il-1,1,-1
          bb(i)=bb(i)/d(i)-e(i)*bb(i+1)
        END DO 
        bb(1:il)=bb(1:il)/SQRT(DOT_PRODUCT(bb(1:il),bb(1:il)))
        IF (ABS(DOT_PRODUCT(aa(1:il),bb(1:il))-1.0_q) < 1.0e-10_q) EXIT
        aa(1:il)=bb(1:il) 
      END DO
! Extract the eigevector of the system from the Lanczos eigenvector
      w=0.0_q   
      DO i=1,il
        w=w+bb(i)*PP(:,:,i)
      END DO
      w=w/SQRT(SUM(w*w))

    RETURN
    END SUBROUTINE eigenvector

!--------------------------------------------------------------------------------------!

    SUBROUTINE eigenvalue(d,e,n)
      INTEGER,INTENT(IN) :: n
      REAL(q),INTENT(INOUT),DIMENSION(n) :: d,e

      INTEGER :: i,itr,l,m
      REAL(q) :: b,c,dd,f,g,p,r,s
      LOGICAL :: small

! Standard QL-factorization

      e(n)=0._q
      DO l=1,n
        itr=0
        DO
          small=.false.
          DO m=l,n-1
            dd=ABS(d(m))+ABS(d(m+1))
            IF (ABS(e(m))+dd == dd) THEN
              small=.true.
              EXIT
            END IF
          END DO
          IF (.not.small) m=n
          IF (m /= l) THEN
            IF (itr == 30) THEN
              IF (iu6 >= 0) WRITE(*,*) 'too many iterations in frequency'
              STOP
            END IF
            itr=itr+1
            g=0.5_q*(d(l+1)-d(l))/e(l)
            r=scaling(g,1._q)
            g=d(m)-d(l)+e(l)/(g+SIGN(r,g))
            s=1._q
            c=1._q
            p=0._q
            DO i=m-1,l,-1
              f=s*e(i)
              b=c*e(i)
              r=scaling(f,g)
              e(i+1)=r
              IF (r == 0._q) THEN
                d(i+1)=d(i+1)-p
                e(m)=0._q  
                EXIT
              END IF
              s=f/r
              c=g/r
              g=d(i+1)-p
              r=(d(i)-g)*s+2._q*c*b
              p=s*r
              d(i+1)=g+p
              g=c*r-b
            END DO
            IF (r==0._q .and. i >= l) CYCLE
            d(l)=d(l)-p
            e(l)=g   
            e(m)=0._q
          END IF
          IF (m == l) EXIT
        END DO
      END DO          
    RETURN
    END SUBROUTINE eigenvalue

!--------------------------------------------------------------------------------------!

    FUNCTION scaling(a,b)
      REAL(q),INTENT(IN) :: a,b
      REAL(q) :: scaling
            
      REAL(q) :: absa,absb,num
            
      absa=ABS(a) 
      absb=ABS(b)
      IF (absa > absb) THEN
        num=absb/absa
        scaling=absa*SQRT(1._q+num*num)
      ELSE
        IF (absb == 0._q) THEN
        scaling=0._q
        ELSE
          num=absa/absb
          scaling=absb*SQRT(1._q+num*num)
        END IF
      END IF
    RETURN      
    END FUNCTION scaling
            
!--------------------------------------------------------------------------------------!

    SUBROUTINE sort(n,ra,rb)
      INTEGER,INTENT(IN) :: n
      REAL(q),INTENT(INOUT),DIMENSION(n) :: ra,rb
            
      INTEGER :: i,ir,j,l
      REAL(q) :: rra,rrb
    
      IF (n < 2) RETURN
      l=n/2+1
      ir=n
      DO
        IF (l > 1) THEN
          l=l-1
          rra=ra(l)
          rrb=rb(l)
        ELSE  
          rra=ra(ir)
          ra(ir)=ra(1)
          rrb=rb(ir)
          rb(ir)=rb(1)  
          ir=ir-1
          IF (ir == 1) THEN
            ra(1)=rra
            rb(1)=rrb
            RETURN
          END IF
        END IF
        i=l
        j=l+l
        DO
          IF (j <= ir) THEN
            IF (j < ir) THEN
              IF (ra(j) < ra(j+1)) j=j+1
            END IF
            IF (rra < ra(j)) THEN
              ra(i)=ra(j)
              rb(i)=rb(j)
              i=j  
              j=j+j
            ELSE
              j=ir+1  
            END IF  
          ELSE
            EXIT 
          END IF
        END DO
        ra(i)=rra
        rb(i)=rrb 
      END DO
    RETURN
    END SUBROUTINE sort

!--------------------------------------------------------------------------------------!

  END MODULE lanczos


