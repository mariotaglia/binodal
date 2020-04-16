subroutine fkfun(x,f,ier)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! User provided routine for kinsol
! x is the input vector
! f is the output vector, kinsol will change x in order to get f = 0
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use results
!use chainsdat
use solver
use system
!use molecules
!use bulk
use const
use results
implicit none
 
integer ntot
real*8 x(3),f(3)
real*16 protemp, protemp1
integer i,j, k, ix, iy, iz, ii, ax, ay, az, temp, iiZ
!real*16 xpotA(dimz),fdisbc,fdisAC
!real*16 xpotB(dimz)
!real*16 psi2(0:dimz+1) ! psi plus boundaries at z=0 and dimz+1
!real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solvent
real*16 m, eta, Penality

! Kinsol
integer*4 ier
real*16 vectalpha(2),vectbeta(2)
real*16 x2alpha,x3alpha,x2beta,x3beta
real*16 mu2alpha,mu2beta,mu3alpha,mu3beta,fealpha,febeta
real*16 potquim2,elib,potquim3
real*16 xsolventalpha,xsolventbeta,xHplusalpha,xHplusbeta
real*16 xOHminalpha,xOHminbeta


! Recovers xh and psi from x
real*16 fracalpha(6), fracbeta(6)

x3alpha=x(1)
x2beta=x(2)
x3beta=x(3)
x2alpha=x2alphafixed

vectalpha(1)=x2alpha
vectalpha(2)=x3alpha
vectbeta(1)=x2beta
vectbeta(2)=x3beta


xsolventalpha=(1.0 -x2alpha-x3alpha)/(1.+KaHplus+KaOHmin)
xsolventbeta=(1.0 -x2alpha-x3alpha)/(1.+KaHplus+KaOHmin)
xHplusalpha=xsolventalpha*KaHplus
xOHminalpha=xsolventalpha*KaOHmin
xHplusbeta=xsolventbeta*KaHplus
xOHminbeta=xsolventbeta*KaOHmin

!print* , xsolventalpha,xsolventbeta
!print* , xHplusalpha,xOHminalpha
!print* , xHplusbeta,xOHminbeta
!print* , x2alpha,x3alpha
!print* , x2beta,x3beta
!stop

call fracasos(vectalpha,fracalpha)
call fracasos(vectbeta,fracbeta)
 
!print* , fracalpha
!print* , fracbeta
!stop 
! Pot quimico respecto de phiA 
Penality=((x2alpha-x2beta)**2+(x3alpha-x3beta)**2)

call mu2(vectalpha,potquim2)
mu2alpha = potquim2

call mu2(vectbeta,potquim2)
mu2beta=potquim2

call mu3(vectalpha,potquim3)
mu3alpha = potquim3

call mu3(vectbeta,potquim3)
mu3beta=potquim3
 
call fe(vectalpha,elib)
fealpha=elib

call fe(vectbeta,elib)
febeta=elib

!print*, mu2beta,mu2alpha,mu3alpha,mu3beta
! ### EQUATIONS TO SOLVE

 f(1)= mu2alpha-mu2beta
 f(1)=f(1)/Penality

 f(2)= mu3alpha-mu3beta
 f(2)= f(2)/Penality

! Recta tangente

 f(3)= (fealpha-febeta&
-((x2alpha-x2beta)/(Ma*vp))*(mu2alpha+mu2beta )/2.&
-((x3alpha-x3beta)/(Mb*vp))*(mu3beta+mu3alpha)/2.) /Penality !	

!print*, f
!stop
!call fracasos(vectalpha,fracalpha)
!call fracasos(vectbeta,fracbeta)
!f(3) = x2alpha/(Ma*vp)+x3alpha/(Mb*vp)+(1-x2alpha-x3alpha)/vs
!f(3) = f(3) -(1/vs)*(log(1-x2alpha-x3alpha))-x2alpha/vp*fracalpha(1)
!f(3) = f(3) - x2beta/(Ma*vp)-x3beta/(Mb*vp)-(1-x2beta-x3beta)/vs
!f(3) = f(3) +(1/vs)*(log(1-x2beta-x3beta))+x2beta/vp*fracbeta(1)
!f(3) = f(3) / Penality

!do iz=1,dimz
!	f(iz+ntot+ntot)= (fealpha-febeta&
!-((x2alpha*mu2beta-x2beta*mu2alpha)/(Ma*vp))&
!-((x3alpha*mu3beta-x3beta*mu3alpha)/(Mb*vp))) /Penality !	
!enddo

iter = iter + 1
norma = 0.0

do i = 1, 3
norma = norma +(f(i))**2    
enddo

ier = 0.0
return
end subroutine
