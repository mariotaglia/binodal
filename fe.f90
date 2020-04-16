subroutine fe(x,elib)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! this routine calculates the free energy of the system and chemical potential 
! of chains
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!use results
!use chainsdat
use system
!use kai
!use bulk
!use molecules
use const
implicit none

integer cc, ccc
real*16 x(2)
real*16 elib,xmsolv,xmOHmin,xmHplus
real*16 Free_energy 
integer i, iz, iiz
!real*8 xtotal(-Xulimit:dimz+Xulimit) ! xtotal for poor solventireal
real*16 xphiA,xphiB,xsolv,fa_A,fa_B,xmphiA,xmphiB,fnc_a,fnc_B,fc_A,fc_B
real*16 frac(6)

! Calculation of xsolvent
xphiA=x(1)
xphiB=x(2)
xsolv=(1.0-xphiA-xphiB)/(1.+KaHplus+KaOHmin)
xmsolv=xsolv/vs
xmHplus=xmsolv*KaHplus
xmOHmin=xmsolv*KaOHmin

 call fracasos(x,frac)
fa_A=frac(1)
fa_B=frac(2)
fnc_A=frac(3)
fnc_B=frac(4)
fc_A=frac(5)
fc_B=frac(6)

xmphiA=xphiA/(Ma*vp)
xmphiB=xphiB/(Mb*vp)

Free_Energy = 0.0

! 1. solvent entropy

  Free_Energy=Free_Energy-xmphiA-xmphiB-xmsolv-xmHplus-xmOHmin

  Free_Energy=Free_Energy +xmphiA*Ma*fa_A+(log(xsolv))/vs

elib=Free_Energy
return 
end subroutine

