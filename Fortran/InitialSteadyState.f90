SUBROUTINE InitialSteadyState

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE
INTEGER		:: iflag,iy
REAL(8)		:: lrhoL,lrhoU
REAL(8), EXTERNAL	:: FnDiscountRate

IF (Display>=1) write(*,*) "Solving for initial steady state"
initialSS = .true.

IF(EquilibriumR==0 .and. CalibrateDiscountRate==0) THEN
	CALL Grids
	CALL IterateBellman
	CALL StationaryDistribution
	
ELSE IF(EquilibriumR==1 .and. CalibrateDiscountRate==0) THEN
	CALL SolveSteadyStateEqum

ELSE IF(CalibrateDiscountRate==1) THEN
	converged = .false.
	neqmiter = 1
	OPEN(3, FILE = trim(OutputDir) // 'DiscountRateCalibration.txt', STATUS = 'replace'); CLOSE(3)
	IF(ngpy==1 .and. deathrate==0.0) THEN !RA model no death
		lrhoL = invlogistic(exp(-0.02_8))
		lrhoU = invlogistic(exp(-0.015_8))
		CALL rtbis(FnDiscountRate,lrhoL,lrhoU,1.0e-6_8,tolrho,rho)
	ELSE
! 		lrhoL = invlogistic(exp(-0.02_8))
! 		lrhoU = invlogistic(exp(-0.01_8))
		lrhoL = invlogistic(exp(-0.016_8))
		lrhoU = invlogistic(exp(-0.014_8))
		CALL rtflsp(FnDiscountRate,lrhoL,lrhoU,1.0e-8_8,tolrho,iflag,maxiterrho)		
	END IF
	converged = .true.	
	
	IF(EquilibriumR==1) CALL SolveSteadyStateEqum
	
END IF

IF(ComputeCumulativeMPC==1) CALL CumulativeConsumption
IF (ComputeDiscountedMPC==1) THEN
	CALL DiscountedMPC
	solnINITSS%mpc = mpc
	solnINITSS%subeff1ass = subeff1ass
	solnINITSS%subeff2ass = subeff2ass
	solnINITSS%wealtheff1ass = wealtheff1ass
	solnINITSS%wealtheff2ass = wealtheff2ass
END IF	 
	
DO iy = 1,ngpy
	ALLOCATE(solnINITSS%A(iy)%val(ACSR(iy)%nz))
	ALLOCATE(solnINITSS%A(iy)%row(ACSR(iy)%n+1))
	ALLOCATE(solnINITSS%A(iy)%col(ACSR(iy)%nz))
	ALLOCATE(solnINITSS%B(iy)%val(BCSR(iy)%nz))
	ALLOCATE(solnINITSS%B(iy)%row(BCSR(iy)%n+1))
	ALLOCATE(solnINITSS%B(iy)%col(BCSR(iy)%nz))
	ALLOCATE(solnINITSS%AU(iy)%val(AUCSR(iy)%nz))
	ALLOCATE(solnINITSS%AU(iy)%row(AUCSR(iy)%n+1))
	ALLOCATE(solnINITSS%AU(iy)%col(AUCSR(iy)%nz))
END DO

solnINITSS%A = ACSR
solnINITSS%B = BCSR
solnINITSS%AU = AUCSR
solnINITSS%V = V
solnINITSS%u = u
solnINITSS%c = c
solnINITSS%s = s
solnINITSS%h = h
solnINITSS%d = d
solnINITSS%bdot = bdot  

solnINITSS%gjoint = gjoint
solnINITSS%gamarg = gamarg
solnINITSS%gbmarg = gbmarg
solnINITSS%gvec = gvec
solnINITSS%gmat = gmat

CALL DistributionStatistics
equmINITSS = EquilibriumType(	ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,labtax,&
								borrwedge,rho,mpshock,prefshock,priceadjust,fundlev,elast,gam,fundbond,profit,dividend,divrate,lumptransfer,equity,caputil,deprec,tfpadj,illassetdrop)
statsINITSS = DistributionStatsType(Ea,Eb,Ec,Elabor,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
									EbN,EbP,Eadjcost,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, &
									Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far)

initialSS = .false.

END SUBROUTINE InitialSteadyState