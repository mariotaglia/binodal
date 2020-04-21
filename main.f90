	!###############################################################################
	!     
	!     Simple brush: Standard Molecular Theory Program 
	!    
	!     Calculates a weak polyelectrolyte brush in poor sv conditions 
	!     Calculates free-energy
	!     Solve BINODAL METHOD
	!     MARCH 2020
	!         
	!###############################################################################
!use pks
use system
use const

	implicit none
	integer i, flagcrash
	real*16 pKw
	!print*, 'Program Simple Brush'
	!print*, 'GIT Version: ', _VERSION
yes=0 ! es para  chequear si encuentra o no xalpha, xbeta
flagcrash=1
	call readinput  ! reads input variables from file
!	call allocation ! allocates memory
	Keo=10**(-pKeo)
	KaA=10**(-pKaA)
	KaB=10**(-pKaB)
	pKw=14.0
	pOHbulk=pKw-pHbulk
	cOHminbulk=10**(-pOHbulk)
	cHplusbulk=10**(-pHbulk)
	!KaHplus=10**(-pkaHplus)
	!KaOHmin=10**(-pkaOHmin)
	!Knew=Keo*(KaB**(KaB/(1+KaB)))*(KaA**(KaA/(1+KaA)))
	!Knew=Keo*KaB*KaA/((KaB+xmphi_OH_free)*(KaA+xmphi_H_free))
	call solve(flagcrash)
	!call fe(cc, ccc)         ! calculates and saves free energy to disk
	call salvar(flagcrash)
	print*, 'Save OK',yes
	call endall     ! clean up and terminate
end
subroutine endall
stop
end subroutine
