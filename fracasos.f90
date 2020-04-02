subroutine fracasos(x,frac)
use system
use const
 implicit none

real*16 frac(6),x(2)
real*16 xphiA,xphiB,xmphiA,xmphiB
real*16 eta,M ,etap
xphiA=x(1);
xphiB=x(2);
xmphiA=xphiA/(Ma*vp);
xmphiB=xphiB/(Mb*vp);

eta=xphiA/xphiB;
etap=xphiB/xphiA;
M=(1.0+KaA)*(1.0+KaB)/(xmphiA*Ma*Knew);
frac(1)=0.5*(1.0+etap+M )-( (0.5*(1.0+etap+M))**2-etap)**0.5;!fA_a
frac(2)=frac(1)*eta;!fB_a
frac(3)=KaA*(1-frac(1))/(1+KaA) ;!fAnc
frac(4)=KaB*(1-frac(2))/(1+KaB) ;!fBnc
frac(5)=(1-frac(1))/(1+KaA);!fAc
frac(6)=(1-frac(2))/(1+KaB);!fBc

!print*,'frac',frac
!stop
!fa=y(1);
!fb=y(2);                                                    
end subroutine
