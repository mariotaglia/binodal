subroutine solve(flagcrash)

!use pks 
use system
use results
use solver
implicit none
real*16 x2betaOK,x2alphaOK
real*16 x2beta0,x2alpha0
real*16 x3betaOK,x3alphaOK
real*16 x3beta0,x3alpha0
real*16 criterio
integer nmax
real*16 xiteri
real*16 xiterii
real*16 xiterj
real*16 xiterjj
real*8 x1(3)
real*8 x1g(3)
integer ier, i,newt, j, ii, jj
integer flagcrash
real*16 potquim2,potquim3,mu2alpha,mu2beta
real*16 xatest(2),xbtest(2),segsolalpha(2),segsolbeta(2)
real*16 Mtemp
real*16 phi0
integer flag
integer ngrid
real*16 arrayalphagrid(2,10000), arraybetagrid(2,10000)
integer gridpoints
real*16 xiter

!   phimin=phiminread; !valor minimo del exponente
!   phimax=phimaxread ! valor maximo del exponente

   nmax=npasos     ! npasos por consola i
   ngrid = npasosgrid
   criterio=1E-3!criterio pra la norma
   gridpoints = 0

! #### GRID SEARCH

linearsolver = 2

! ida

   do i=1,ngrid
   do ii=1,ngrid
   do j=1,ngrid
   do jj=1,ngrid
  
      xiteri=phimin+float(i-1)*(phimax-phimin)/float(ngrid-1)
      xiterii=phimin+float(ii-1)*(phimax-phimin)/float(ngrid-1)
      xiterj=phimin+float(j-1)*(phimax-phimin)/float(ngrid-1)
      xiterjj=phimin+float(jj-1)*(phimax-phimin)/float(ngrid-1)
   
      x2alphafixed= 10**xiteri ! x2phialpha
      x3alphaOK = 10**xiterj
      x2betaOK = 10**xiterii
      x3betaOK = 10**xiterjj

      x1(1)=x3alphaOk
      x1g(1)=x1(1)
      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)


!      print*, x2alphafixed, x1(1), x1(2), x1(3)
      call call_kinsol(x1, x1g, ier)

       if (norma.lt.criterio) then !esto es para saber si encontró o no solución 
         gridpoints = gridpoints + 1
         print*,'Grid Point OK',gridpoints
         print*,'x2alpha,x3alpha,x2beta,x3beta',x2alphafixed,x1(1),x1(2),x1(3)

         arrayalphagrid(1,gridpoints)=x2alphafixed
         arrayalphagrid(2,gridpoints)=x1(1)
         arraybetagrid(1,gridpoints)=x1(2)
         arraybetagrid(2,gridpoints)=x1(3)

       endif

   enddo
   enddo
   enddo
   enddo

linearsolver = 2

! ### Seach using grid as seed

conteo=0  !lo uso  para  guardar en el arrayalpha y arraybeta

do ii = 1, gridpoints

phi0 = log10(arrayalphagrid(1,ii))
x3alphaOK = arrayalphagrid(2,ii)
x2betaOK = arraybetagrid(1,ii)
x3betaOK = arraybetagrid(2,ii)

! ida

   do i=1,nmax
  
      xiter=phi0+float(i-1)*(phimax-phi0)/float(nmax-1)
   
      x2alphafixed= 10**xiter ! x2phialpha

      x1(1)=x3alphaOk
      x1g(1)=x1(1)
      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)

      call call_kinsol(x1, x1g, ier)
      flag = 0
       if (norma.lt.criterio) then !esto es para saber si encontró o no solución 
         flag = 1
         conteo=conteo+1
         yes=yes+1
         print*,'Yes',yes
         print*,'x2alpha,x3alpha,x2beta,x3beta',x2alphafixed,x1(1),x1(2),x1(3)

         x3alphaOK= x1(1)
         x2betaOK = x1(2)
         x3betaOK = x1(3)

         arrayalpha(1,conteo)=x2alphafixed
         arrayalpha(2,conteo)=x1(1)
         arraybeta(1,conteo)=x1(2)
         arraybeta(2,conteo)=x1(3)

       endif
   enddo


!!! VUELVE

phi0 = log10(arrayalphagrid(1,ii))
x3alphaOK = arrayalphagrid(2,ii)
x2betaOK = arraybetagrid(1,ii)
x3betaOK = arraybetagrid(2,ii)

   do i=1, nmax

      xiter=phi0+float(i-1)*(phimin-phi0)/float(nmax-1)
      x2alphafixed= 10**xiter ! x2phialpha

      print*, i, x2alphafixed

      x1(1)=x3alphaOk
      x1g(1)=x1(1)
      x1(2)=x2betaOK     !x2phibeta inicial
      x1g(2)=x1(2)
      x1(3)=x3betaOK     !x3phibeta inicial
      x1g(3)=x1(3)

      call call_kinsol(x1, x1g, ier)

      flag = 0

       if (norma.lt.criterio) then !esto es para saber si encontró o no solución 
         flag = 1
         conteo=conteo+1
         yes=yes+1
         print*,'Yes',yes
         print*,'x2alpha,x3alpha,x2beta,x3beta',x2alphafixed,x1(1),x1(2),x1(3)

         x3alphaOK= x1(1)
         x2betaOK = x1(2)
         x3betaOK = x1(3)

         arrayalpha(1,conteo)=x2alphafixed
         arrayalpha(2,conteo)=x1(1)
         arraybeta(1,conteo)=x1(2)
         arraybeta(2,conteo)=x1(3)
       endif
   enddo


enddo ! grid search 

return
end subroutine

