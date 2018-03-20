SUBROUTINE ImpulseResponses

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: it,it1
CHARACTER	:: IRFDir*20

!allocate arrays
CALL AllocateArrays

stickytransition = .false.

!set up deltatransvec
CALL PowerSpacedGrid (Ttransition,deltatransparam,deltatransmin,deltatransmax,deltatransvec)

nendtrans = min(50,Ttransition)
cumdeltatrans(1) = deltatransvec(1)
DO it = 2,Ttransition
	cumdeltatrans(it) = cumdeltatrans(it-1) + deltatransvec(it)
END DO
OPEN(3, FILE = trim(OutputDir) // 'deltatransvec.txt', STATUS = 'replace'); CALL WriteMatrix(3,Ttransition,1,deltatransvec)



!Monetary policy shock
IF(IncludeMonetaryShock==1) THEN
	IF(Display>=1) write(*,*)'Solving for monetary policy shock IRF'	
	irfpointer => irfstruct(0)

	equmTRANS(:) = equmINITSS		
	equmTRANS(1)%mpshock = equmINITSS%mpshock + MonetaryShockSize
	DO it = 2,Ttransition-nendtrans
		equmTRANS(it)%mpshock =equmINITSS%mpshock *(1.0-MonetaryShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%mpshock * (MonetaryShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%mpshock = equmINITSS%mpshock

	IRFDir = "Monetary"
	CALL IRFSequence(IRFDir)

END IF

!Forward Guidance shock
IF(IncludeForwardGuideShock==1) THEN
	forwardguide = .true.
	IF(Display>=1) write(*,*)'Solving for forward guidance shock IRF'	
	irfpointer => irfstruct(0)

	equmTRANS(:) = equmINITSS
	it1 = MINLOC(cumdeltatrans, 1, MASK = cumdeltatrans>=ForwardGuideShockQtrs)
	equmTRANS(1:it1-1)%mpshock = equmINITSS%mpshock
	equmTRANS(it1)%mpshock = equmINITSS%mpshock + ForwardGuideShockSize
	DO it = it1+1,Ttransition-nendtrans
		equmTRANS(it)%mpshock =equmINITSS%mpshock *(1.0-MonetaryShockPers**deltatransvec(it-1)) + equmTRANS(it-1)%mpshock * (MonetaryShockPers**deltatransvec(it-1))

	END DO	
	equmTRANS(Ttransition-nendtrans+1:Ttransition)%mpshock = equmINITSS%mpshock

	IRFDir = "ForwardGuide"
	CALL IRFSequence(IRFDir)
	forwardguide = .false.

END IF


END SUBROUTINE ImpulseResponses