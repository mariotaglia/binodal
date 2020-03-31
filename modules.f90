module system
!real*16 delta   ! delta is the discretization lenght in z direction
!real*16 sigmaA
real*16 Ma
integer yes
real*16 Mb
integer npasos
integer npasosgrid
real*16 phimin, phimax
real*16 sigmaB
real*16 csalt
real*16 nor
!real*16 pKaA
!real*16 pkaB
!real*16 pKaANa
!real*16 pkaBCl
real*16 pkEo
real*16 Keo
!real*16 pHbulk
!real*8 st
!integer VUELTA
endmodule

module results
real*16 x2alphafixed
!real*16 x3alphafixed
real*16 arrayalpha(2,4000000)
real*16 arraybeta(2,4000000)
integer conteo
endmodule

module const
!real*16, parameter :: pi = 3.14159 ! pi 
!real*16, parameter :: Na = 6.02d23 ! Avogadro's number
!real*16, parameter :: lb = 0.714   ! bjerrum lenght in water in nm
real*16, parameter :: vs = 0.03  ! bjerrum lenght in water in nm
real*16, parameter :: vp = 0.11   ! bjerrum lenght in water in nm

!real*16 constq
!real*16 pKw
endmodule


module solver
!real*16, allocatable :: xflag(:)
!integer infile
!real*16, parameter :: error = 1.0d-6 ! maximum kinsol norm
real*16 norma
integer iter
integer linearsolver
endmodule


!module molecules
!real*16 zpos, zneg, zpolA, zpolB ! charges of cation, anions and polyelectrolyte segment
!real*16 vsalt, vpol   ! volume of salt and polyelectrolyte segments in units of vsol
!real*16 vsol             ! solvent volume 
!real*16 K0A, K0B ,K0ANa,K0BCl, K0Eo!K0
!endmodule

!module bulk
!real*16 xHplusbulk, xOHminbulk ! bulk volume fraction of H+ and OH-
!real*16 xposbulk, xnegbulk     ! bulk volume fraction of cation and anion
!real*16 xsolbulk               ! bulk volume fraction of solvent
!real*16 expmupos, expmuneg, expmuHplus, expmuOHmin  ! exp(-beta*mu)*(bulk volume fraction), where mu is the chemical potential
!endmodule

!module rand
!integer seed
!endmodule

module mkinsol
double precision, allocatable :: pp(:)
endmodule
