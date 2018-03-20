SUBROUTINE SetParameters

USE Parameters
USE Globals
USE Procedures

IMPLICIT NONE

INTEGER 	:: iy
REAL 		:: la,lb,lc

!OUTPUT DIR
OutputDir =	"~/FortranOutputDir/BaselineOutputSubdir/"
EarningsProcessDir	= "earnings_input"

CALL system ("mkdir -p " // trim(OutputDir))

!OPTIONS

!display options
Display              		= 1
ReportNonMonotonicity   	= 0

!run options
CalibrateDiscountRate	= 0 
EquilibriumR		 	= 1
ComputeCumulativeMPC 	= 1
DoImpulseResponses 		= 1
DoPriceExperiments		= 1
SaveTime1PolicyFns 		= 1
SaveCumPolicyFnsIRF 	= 0
ComputeDiscountedMPC 	= 1

!labor risk options
ReadEarningsProcess 	= 1
NoLaborSupply		= 0	!only one of the labor supply options
LaborSupplySep		= 1	
LaborSupplyGHH 		= 0
ScaleDisutilityIdio	= 0
ImposeMaxHours 		= 1
PerfectAnnuityMarkets	= 1
GovBondResidualZeroWorld= 1	!imposes closed economy and imputes residual bond holdings to govt
AdjustProdGridFrisch 	= 1
adjfricshgridfrac = 0.85 !fraction of Frisch to adjust by

!calibration options
CalibrateCostFunction		= 0
CalibrateRhoAtInitialGuess  = 0!1
MatchRelativeToTargetOutput	= 0
ImposeEqumInCalibration 	= 0

!adjustment cost options
OneAssetNoCapital		= 0
PinKappa1ByKappa02 		= 1!0

!transition computation options
SolveStickyPriceTransition	= 1
ConvergenceRelToOutput 		= 1
FirmDiscountRate			= 5	!1 for rho, 2 for rb initial steady state, 3 for ra initial steady state, 4 for rb transition, 5 for ra transition
bondelastrelgdp 			= 1.0 !bigger for smaller interest rate movements, closer to zero for larger interest rate movements. relative to steady state gdp
bondadjust 					= 0.1 !more responsive interest rate when closer to zero

!dividend options
DividendFundLumpSum 			= 1 
DistributeProfitsInProportion 	= 1 !distributes profdistfrac of profits as dividends, and (1-profdistfrac) to households in proportion to productivity
profdistfrac 					= 0.33 !set to alpha to neutralize effect of profits on redistribution between liquid and illiquid assets
TaxHHProfitIncome 				= 1 !taxes profit income at labor tax rate if DistributeProfitsInProportion = 1


!government bc options
AdjGovBudgetConstraint 		= 2 !1 for adjust spending, 2 for adjust lump sum taxes, 3 for let debt adjust (choose options below for financing), 4 for adjust proportional tax
GovExpConstantFracOutput 	= 0 !only active if AdjGovBudgetConstraint==3
taxincrstart 		= 1 !quarters after shock that fiscal policy adjusts
taxincrdecay 		= 0.02 !decay rate for tax increase higher for faster decay


!MONETARY POLICY SHOCK
IncludeMonetaryShock	= 1
MonetaryShockSize 		= -0.0025 !percentage points
MonetaryShockPers 		= exp(-0.5) !0.5 !quarterly

!FORWARD GUIDANCE SHOCK
IncludeForwardGuideShock= 0!1
ForwardGuideShockSize 	= -0.0025 !percentage points
ForwardGuideShockPers 	= exp(-0.5) !quarterly
ForwardGuideShockQtrs 	= 9 !number of quarters in advance (set phifg below)


!CALIBRATION OPTIONS
EstimateKappa0		= 1
EstimateKappa1		= 0
EstimateKappa2		= 1
EstimateKappa3		= 0
EstimateKappa4		= 0
EstimateRho			= 1
EstimateBorrWedge	= 1
EstimateGamma		= 0


MatchMeanIll		= 0
MatchKYratio		= 1
MatchMedianIll		= 0
MatchP75Ill			= 0
MatchFracIll0		= 0
MatchMeanLiq		= 1
MatchMedianLiq		= 0
MatchFracLiq0		= 1
MatchFracLiqNeg		= 1
MatchFracIll0Liq0 	= 1!0

defnbclose			= 0.0
defnaclose 			= 0.0

ndfls		= 2

!SOLUTION PARAMETERS
maxiter 		= 500
Vtol			= 1.0e-8
maxiterKFE		= 2000
KFEtol			= 1.0e-12
deltass  		= 1.0e6
deltakfe 		= 1.0e6
dVamin 			= 1.0e-8
dVbmin 			= 1.0e-8

tolequmss		= 1.0e-10
stepequmss		= 0.10
maxiterequmss	= 40
maxiterrho 		= 30
tolrho			= 1.0e-8

toltransition	= 1.0e-6
deltatransmin	= 1.0/3.0
deltatransmax	= 40.0
deltatransparam	= 0.35
maxitertranssticky	= 5000
stepstickytransK  = 0.01
stepstickytransB  = 0.001 !0.001		

deltacumcon = 0.01 !deltatransmin !0.01 !set to a low number like 0.01 for accurate steady state MPCs, and to deltratransmin for IRF consistency

!discount rates
rho		=  0.01272892513 !0.01272892513(baseline), 0.0133348071(frisch=0.5), 0.0123887061(frisch=1.5), 0.0125456234 (ghh)


!preferences
deathrate	= 1.0/(4.0*45.0) !poisson death rate
gam			= 1.0	!risk aversion
prefshock	= 1.0


!liquid assets
rb			= 0.02/4.0 !liquid return
borrwedge 	= 0.0148846 !0.019663 ! !quarterly wedge between rb and rborr: intermediation cost  
borrwedgemax= 0.09
blim 		 = -1.0 	!borrowing limit multiple of quarterly output

rborr = rb + borrwedge


!withdrawal costs
kappa0_w	= 0.04383
kappa2_w	= 0.40176
kappa3		= 0.03 *2.92/4.0 !as a fraction of average illiquid assets
kappa4_w 	= 0.0

IF(PinKappa1ByKappa02==0) kappa1_w	=  0.48236
IF(PinKappa1ByKappa02==1) kappa1_w	= ((1.0-kappa0_w)*(1.0+kappa2_w))**(-1.0/kappa2_w)

kappa2min   = 0.05 !to make sure there is enough curvature for calibration

!deposit costs
kappa0_d = kappa0_w
kappa1_d = kappa1_w
kappa2_d = kappa2_w
kappa4_d = kappa4_w
	

theta 		= 100.0 !price adjustment cost
phitaylor 	= 1.25  ! inflation coefficient in taylor rule, 0.0 for fixed nominal interest rate at ss.
phifg 		= phitaylor !with forward guidance use phifg instead of phitaylor pre-shock, and phitaylor after (not equal 1)


pi 		= 0.0		!inflation rate: this is steady state target	
rnom = rb + pi		!nominal interest rate (fisher equation): this will be constant in taylor rule
mpshock 	= 0.0		!shock to nominal interest rate in taylor rule

tfp 	= 1.0
elast 	= 10.0 !elasticity of DS aggregator
gap 	= 0.0 !steady state output gap
mc = (elast-1.0)/elast

alpha 		= 0.33
alphatilde 	= (alpha**alpha) * ((1.0-alpha)**(1.0-alpha))
deprec 		= 0.07/4.0	!depreciation rate

fundlev 	= 0.0
utilelast 	= 0.0
utilelastalpha  = 1.0 + utilelast-alpha*utilelast

!government
labtax 		= 0.30
lumptransferpc = 40000*labtax/(115000.0*2.92*4.0)
! lumptransfer = 0.10
corptax 		= 0.0
ssdebttogdp 	= 0.26*4 !if foreign sector assumed to hold residual bonds

!calibration targets (relative to quarterly output)
targetMeanIll 		= 2.92 * 4.0
targetMeanLiq  		= 0.2* 4.0
targetMedianIll 	= 0.21 * targetMeanIll
targetP75Ill 		= 0.71 * targetMeanIll
targetMedianLiq 	= 0.085 * targetMeanLiq
targetFracIll0 		= 0.115 + 0.12
targetFracLiq0 		= 0.35
targetFracLiqNEG	= 0.15
targetFracIll0Liq0 	= 0.10

lumptransfer = lumptransferpc*targetMeanIll


IF (DividendFundLumpSum ==0) THEN
	targetKYratio 	= targetMeanIll/(1.0 - fundlev)

ELSE IF (DividendFundLumpSum ==1) THEN
	
	la = -(deprec+rb*fundlev) 
	IF(DistributeProfitsInProportion==0) lb = ((elast-1.0)/elast)*alpha + (1.0/(1.0-fundlev))*targetMeanIll*(deprec+rb*fundlev) + (1.0/elast)*(1.0-corptax)
	IF(DistributeProfitsInProportion==1) lb = ((elast-1.0)/elast)*alpha + (1.0/(1.0-fundlev))*targetMeanIll*(deprec+rb*fundlev) + (1.0/elast)*(1.0-corptax)*profdistfrac
	lc = -(1.0/(1.0-fundlev)) * targetMeanIll *((elast-1)/elast) *alpha
	
	targetKYratio = (-lb+sqrt(lb**2-4*la*lc)) / (2*la)
END IF

!if solving for equilibrium, these are guesses
KYratio = targetKYratio 
KNratio = (tfp*KYratio)**(1.0/(1.0-alpha))
rcapital = mc*alpha/KYratio
wage = mc*(1.0-alpha)*tfp*(KNratio**alpha)
netwage = (1.0-labtax)*wage
IF(DividendFundLumpSum==1) divrate = 0.0
IF(DividendFundLumpSum==0) divrate =  (1.0-corptax)*(1.0-mc)/KYratio !outside of steady state include price adjustments
IF(DistributeProfitsInProportion==1) divrate =  profdistfrac*divrate

ra = (rcapital - deprec + divrate - fundlev*rb)/(1.0-fundlev)

IF(NoLaborSupply==1) THEN
	!scale efficiency units so that output euqals 1
	meanlabeff = KYratio**(alpha/(alpha-1))
ELSE
	!scale efficiency units so that average hours would be 1/3 in eqm if cov(z,h)=0
	meanlabeff = (KYratio/KNratio)/(1.0/3.0)
END IF		

frisch 		= 1.0 	!frisch elasticity labor supply

IF(OneAssetNoCapital==1) THEN 
	alpha = 0.0
	profdistfrac = 0.0
	rcapital = 0.0
	wage = mc*(1.0-alpha)*tfp
	netwage = (1.0-labtax)*wage
	ra = 0.0
	meanlabeff = 3.0
	kappa0_w = 100.0
	kappa0_d = 100.0
	FirmDiscountRate = 4
	
END IF

!guess chi's so that at average wages and average consumption hours =1/3 (sets C/Y = 0.75)
IF (NoLaborSupply==1)	chi	= 0.0 
IF (LaborSupplySep==1)	chi	= meanlabeff / (0.75 **(-gam) * (1.0/3.0)**(1.0/frisch))
IF (LaborSupplyGHH==1)	chi = meanlabeff / ((1.0/3.0)**(1.0/frisch)) 

!requires an intial guess of profits if DistributeProfitsInProportion = 1, since output not known in advance
IF(DistributeProfitsInProportion==1) THEN
	!these will be updated in iterations
	IF(OneAssetNoCapital==0) THEN 
		IF (LaborSupplySep==1) profit = (1.0-mc)*16.0/KYratio
		IF (LaborSupplyGHH==1) profit = (1.0-mc)*14.0/KYratio
		IF (NoLaborSupply==1) profit = (1.0-mc)*12.0/KYratio
	ELSE IF(OneAssetNoCapital==1) THEN 
		profit = (1.0-mc)*1.0
	END IF

END IF	



!allocate large arrays
ALLOCATE(V(ngpa,ngpb,ngpy),Vnew(ngpa,ngpb,ngpy))
ALLOCATE(u(ngpa,ngpb,ngpy),c(ngpa,ngpb,ngpy),h(ngpa,ngpb,ngpy),d(ngpa,ngpb,ngpy),s(ngpa,ngpb,ngpy),bdot(ngpa,ngpb,ngpy))
ALLOCATE(ccum1(ngpa,ngpb,ngpy),ccum2(ngpa,ngpb,ngpy),ccum4(ngpa,ngpb,ngpy))
ALLOCATE(dcum1(ngpa,ngpb,ngpy),dcum2(ngpa,ngpb,ngpy),dcum4(ngpa,ngpb,ngpy))
ALLOCATE(gjoint(ngpa,ngpb,ngpy),gamarg(ngpa,ngpy),gbmarg(ngpb,ngpy),gvec(naby),gmat(nab,ngpy),gabmarg(ngpa,ngpb),gabcum(ngpa,ngpb))
ALLOCATE(mpc(ngpa,ngpb,ngpy),subeff1ass(ngpa,ngpb,ngpy),subeff2ass(ngpa,ngpb,ngpy),wealtheff1ass(ngpa,ngpb,ngpy),wealtheff2ass(ngpa,ngpb,ngpy))

!allocate solution types
CALL AllocateSolutionType(solnINITSS)
CALL AllocateSolutionType(solnFINALSS)

!allocate cumulative policy types
CALL AllocateCumulativePolicyType(cumINITSS)

!allocate COO matrices
DO iy = 1,ngpy
	ALLOCATE(ACOO(iy)%val(5*nab),ACOO(iy)%row(5*nab),ACOO(iy)%col(5*nab))
	ALLOCATE(AUCOO(iy)%val(3*nab),AUCOO(iy)%row(3*nab),AUCOO(iy)%col(3*nab))
END DO





calibrating = .false.
iteratingtransition = .false.
forwardguide = .false.

END SUBROUTINE SetParameters
