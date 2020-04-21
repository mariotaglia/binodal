
subroutine mu3(x,potquim3)
use const
use system
implicit none
real*16 x(2)
real*16 potquim3,fa_A,fa_B,fc_a,fc_b,fnc_A,fnc_b
real*16 xphiA,xphiB,xmphiA,xmphiB,xsolv,xOHmin
real*16 frac(6)

xphiA=x(1)
xphiB=x(2)
xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)
xsolv=(1.0-xphiA-xphiB)/(1.+cHplusbulk+cOHminbulk)
xOHmin=xsolv*cOHminbulk
call fracasos(x,frac)
fa_A=frac(1)
fa_B=frac(2)
fnc_a=frac(3)
fnc_b=frac(4)
fc_A=frac(5)
fc_B=frac(6)
potquim3= log(xmphiB*vs)-(Mb*vp/vs)*log(xsolv)+Mb*(log(fc_B)-log(xOHmin)+log(xsolv) )
!return (potquim3)
!print* ,potquim3
!stop
end subroutine 
