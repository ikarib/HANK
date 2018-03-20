SUBROUTINE Calibration

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER					:: ip,npt,maxfun,iprint,j,iflag
REAL(8)					:: lrhobeg, lrhoend,lrhoL,lrhoU
REAL(8), ALLOCATABLE	:: x(:),w(:)
REAL(8), EXTERNAL 		:: FnDiscountRate

calibrating = .true.

IF(CalibrateRhoAtInitialGuess==1) THEN
	IF (Display>=1) write(*,*) "Calibrating rho at initial steady state"
	converged = .false.
	neqmiter = 1
	OPEN(3, FILE = trim(OutputDir) // 'DiscountRateAtInitialGuess.txt', STATUS = 'replace'); CLOSE(3)
	lrhoL = invlogistic(exp(-0.02_8))
	lrhoU = invlogistic(exp(-0.01_8))
	CALL rtflsp(FnDiscountRate,lrhoL,lrhoU,1.0e-8_8,tolrho,iflag,maxiterrho)
	converged = .true.
END IF



!set up parameters
nparam = 0
IF (EstimateKappa0 == 1) 		nparam = nparam+1
IF (EstimateKappa1 == 1) 		nparam = nparam+1
IF (EstimateKappa2 == 1) 		nparam = nparam+1
IF (EstimateKappa3 == 1) 		nparam = nparam+1
IF (EstimateKappa4 == 1) 		nparam = nparam+1
IF (EstimateRho == 1) 			nparam = nparam+1
IF (EstimateBorrWedge == 1) 	nparam = nparam+1
IF (EstimateGamma == 1) 	nparam = nparam+1

IF (ALLOCATED(paramguess)) DEALLOCATE(paramguess)
IF (ALLOCATED(paramout)) DEALLOCATE(paramout)
IF (ALLOCATED(paramscale)) DEALLOCATE(paramscale)
IF (ALLOCATED(paramub)) DEALLOCATE(paramub)
IF (ALLOCATED(paramlb)) DEALLOCATE(paramlb)
IF (ALLOCATED(x)) DEALLOCATE(x)
ALLOCATE(paramguess(nparam),paramout(nparam),paramscale(nparam),paramub(nparam),paramlb(nparam),x(nparam))

!set up moments
nmoments = 0
IF (MatchMeanIll == 1)		nmoments = nmoments+1
IF (MatchKYratio == 1)		nmoments = nmoments+1
IF (MatchMedianIll == 1)	nmoments = nmoments+1
IF (MatchP75Ill == 1)		nmoments = nmoments+1
IF (MatchFracIll0 == 1)		nmoments = nmoments+1
IF (MatchMeanLiq == 1)		nmoments = nmoments+1
IF (MatchMedianLiq == 1)	nmoments = nmoments+1
IF (MatchFracLiq0 == 1)		nmoments = nmoments+1
IF (MatchFracLiqNeg == 1)	nmoments = nmoments+1
IF (MatchFracIll0Liq0 == 1)	nmoments = nmoments+1

!set up weighting matrix and n vector
IF (ALLOCATED(diagweight)) DEALLOCATE(diagweight)
IF (ALLOCATED(nobsvec)) DEALLOCATE(nobsvec)
ALLOCATE(diagweight(nmoments),nobsvec(nmoments))
diagweight = 1.0

!make guess: with bound constraints, scale is set in one go below, based on bounds. the initial starting points are the midpoints of the bounds
ip = 0
IF(EstimateKappa0==1) THEN
    ip=ip+1
	paramguess(ip) = invlogistic(kappa0_w)
	paramub(ip) = 0.10 !0.95
	paramscale(ip) = 1.0/(invlogistic(paramub(ip))-paramguess(ip))
END IF
IF(EstimateKappa1==1) THEN
    ip=ip+1
	paramguess(ip) = log(kappa1_w)
	paramub(ip) = 3.0
	paramscale(ip) = 1.0/(log(paramub(ip))-paramguess(ip))
END IF
IF(EstimateKappa2==1) THEN
    ip=ip+1
	paramguess(ip) = log(kappa2_w-kappa2min)
	paramub(ip) = 4.0
	paramscale(ip) = 1.0/(log(paramub(ip)-kappa2min)-paramguess(ip))
END IF
IF(EstimateKappa3==1) THEN
    ip=ip+1
	paramguess(ip) = log(kappa3)
	paramub(ip) = 10.0
	paramscale(ip) = 1.0/(log(paramub(ip))-paramguess(ip))
END IF
IF(EstimateKappa4==1) THEN
    ip=ip+1
	paramguess(ip) = log(kappa4_w)
	paramub(ip) = 0.5
	paramscale(ip) = 1.0/(log(paramub(ip))-paramguess(ip))
END IF
IF(EstimateRho==1) THEN
    ip=ip+1
	paramguess(ip) = invlogistic(exp(-rho))
	paramub(ip) = 0.02
	paramscale(ip) = 1.0/(invlogistic(exp(-paramub(ip)))-paramguess(ip))
END IF
IF(EstimateBorrWedge==1) THEN
    ip=ip+1
	paramguess(ip) = invlogistic(borrwedge/borrwedgemax)
	paramub(ip) = 0.95*borrwedgemax
	paramscale(ip) = 1.0/(invlogistic(paramub(ip)/borrwedgemax)-paramguess(ip))

END IF
IF(EstimateGamma==1) THEN
    ip=ip+1
	paramguess(ip) = log(gam)
	paramub(ip) = 3.0
	paramscale(ip) = 1.0/(log(paramub(ip))-paramguess(ip))
END IF


!write paramguess to file
OPEN(3, FILE = trim(OutputDir) // 'paramguess' //   '.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,Nparam,1,paramguess)
OPEN(3, FILE = trim(OutputDir) // 'paramlb' //   '.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,Nparam,1,paramlb)
OPEN(3, FILE = trim(OutputDir) // 'paramub' //   '.txt', STATUS = 'replace'); CALL WriteMatrixLong(3,Nparam,1,paramub)



x = paramguess*paramscale
npt = 2*nparam+1
lrhobeg = 1.0
lrhoend = 1.0e-5
maxfun = 500*(nparam+1) 
iprint = 3
IF(ALLOCATED(w)) DEALLOCATE(w)
ALLOCATE(w((npt+11)*(npt+nparam) +nparam*(3*nparam+11)/2) )

!do DFLS/DFBOLS minimization: objective function in dfovec.f90
!NB: DFBOLS automatically makes points the midpoints of the bounds
!DFLS uses the actual guesses with no bounds
OPEN(4, FILE = trim(OutputDir) // 'iterations' //   '.txt', STATUS = 'replace')
DO j = 1,ndfls
	
	IF(CalibrateRhoAtInitialGuess==1) THEN
		IF (Display>=1) write(*,*) "Calibrating rho at initial steady state"
		converged = .false.
		neqmiter = 1
		OPEN(3, FILE = trim(OutputDir) // 'DiscountRateAtInitialGuess.txt', STATUS = 'replace'); CLOSE(3)
		lrhoL = invlogistic(exp(-0.02_8))
		lrhoU = invlogistic(exp(-0.01_8))
		CALL rtflsp(FnDiscountRate,lrhoL,lrhoU,1.0e-8_8,tolrho,iflag,maxiterrho)
		converged = .true.
	END IF
	
	write(4,*) '********************************** '
	write(4,*) 'DFLS MINIMIZATION ATTEMPT ', j
	CALL NEWUOA_H(nparam,npt,x,lrhobeg,lrhoend,iprint,maxfun,w,nmoments)
	write(4,*) '********************************** '
END DO
CLOSE(4)
paramout = x(1:nparam)/paramscale


!extract parameters at solution
ip = 0
IF(EstimateKappa0==1) THEN
    ip=ip+1
	kappa0_w = logistic(paramout(ip))
	write(*,*) ' kappa0_w sol: ',kappa0_w
END IF

IF(EstimateKappa1==1) THEN
    ip=ip+1
    kappa1_w = exp(paramout(ip))
	write(*,*) ' kappa1_w sol: ',kappa1_w
END IF

IF(EstimateKappa2==1) THEN
    ip=ip+1
    kappa2_w = kappa2min + exp(paramout(ip))
	write(*,*) ' kappa2_w sol: ',kappa2_w
END IF

IF(EstimateKappa3==1) THEN
    ip=ip+1
    kappa3 = exp(paramout(ip))
	write(*,*) ' kappa3 sol: ',kappa3
END IF

IF(EstimateKappa4==1) THEN
    ip=ip+1
    kappa4_w = exp(paramout(ip))
	write(*,*) ' kappa4_w sol: ',kappa4_w
END IF

IF(EstimateRho==1) THEN
    ip=ip+1
	rho = -log(logistic(paramout(ip)))
	write(*,*) ' rho sol: ',rho
END IF

IF(EstimateBorrWedge==1) THEN
    ip=ip+1
	borrwedge = borrwedgemax*logistic(paramout(ip))	
	write(*,*) ' borrwedge sol: ',borrwedge
END IF

IF(EstimateGamma==1) THEN
    ip=ip+1
    gam = exp(paramout(ip))
	write(*,*) ' gam sol: ',gam
END IF

IF(PinKappa1ByKappa02==1)THEN
	kappa1_w = ((1.0-kappa0_w)*(1.0+kappa2_w))**(-1.0/kappa2_w)
END IF	

kappa0_d = kappa0_w
kappa1_d = kappa1_w
kappa2_d = kappa2_w
kappa4_d = kappa4_w


calibrating = .false.

!implied aggregate statistics: note that lump transfer is based on output=1.5, not actual output
bond = Eb
investment = deprec*capital
priceadjust = 0.0
profit = (1.0-mc)*capital/KYratio - priceadjust
IF(DistributeProfitsInProportion==0) dividend = profit*(1.0-corptax)
IF(DistributeProfitsInProportion==1) dividend = profdistfrac*profit*(1.0-corptax)
output = tfp*(capital**alpha)*(labor**(1.0-alpha))
fundbond = -capital*fundlev
bondelast = bondelastrelgdp*output
caputil 	= 1.0
tfpadj = ((tfp**(1.0+utilelast)) * (mc*alpha/rcapital)**(alpha*utilelast))**(1.0/utilelastalpha)
taxrev = labtax*wage*labor - lumptransfer + corptax*profit
IF(DistributeProfitsInProportion == 1 .and. TaxHHProfitIncome == 1) taxrev = taxrev + labtax*(1.0-profdistfrac)*profit*(1.0-corptax)


IF(GovBondResidualZeroWorld==0) THEN
	govbond = -ssdebttogdp*output
	govexp = taxrev + rb*govbond 
	worldbond = -bond-govbond-fundbond
ELSE IF(GovBondResidualZeroWorld==1) THEN
	worldbond = 0.0
	govbond = -bond-worldbond-fundbond
	govexp = taxrev + rb*govbond		
END IF

END SUBROUTINE Calibration
