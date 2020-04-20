subroutine fracasos(x,frac)
use system
use const
 implicit none
integer i
real*16 frac(6),x(2)
real*16 xphiA,xphiB,xmphiA,xmphiB
real*16 eta,M ,etap,Kap,Kbp,Kop
xphiA=x(1);
xphiB=x(2);
xmphiA=xphiA/(Ma*vp);
xmphiB=xphiB/(Mb*vp);

!print* ,x
eta=xphiB/xphiA;
etap=xphiA/xphiB;
Kap=KaA/KaHplus
Kbp=KaB/KaOHmin
Kop=Keo*xmphiA*Ma*vs

Knew=Kop*Kap*Kbp/((1.0+Kap)*(1.0+Kbp))
frac(2)=0.5*(1.0+etap+etap/Knew )-( (0.5*(1.0+etap+etap/Knew))**2-etap)**0.5;!fB_a
frac(1)=frac(2)*eta;!fA_a
frac(3)=(1-frac(1))/(1+Kap) ;!fAnc
frac(4)=(1-frac(2))/(1+Kbp) ;!fBnc
frac(5)=(1-frac(1))*Kap/(1+Kap);!fAc
frac(6)=(1-frac(2))*Kbp/(1+Kbp);!fBc

!print*, eta,etap
!print* ,frac
!stop
!do i = 1,6
!  if (frac(i) .lt. 10**(-15))then
!    frac(i)=10**(-15)
! endif
!enddo
                                  
end subroutine
