MODULE Globals
USE Parameters

IMPLICIT NONE

!GLOBALS FOR DIRECTORIES
character(len=100)  OutputDir,OutputDirIRF,EarningsProcessDir
character(len=100)  InputParamFile

!OPTIONS GLOBALS
integer	:: Display, ReportNonMonotonicity, NoLaborSupply,LaborSupplySep,LaborSupplyGHH,EquilibriumR, CalibrateCostFunction, DoImpulseResponses,CalibrateDiscountRate,SaveCumPolicyFnsIRF,ComputeDiscountedMPC
integer	:: SolveStickyPriceTransition, FirmDiscountRate, ImposeEqumInCalibration
integer :: ComputeCumulativeMPC, ScaleDisutilityIdio,SaveTime1PolicyFns,DoPriceExperiments,DistributeProfitsInProportion,TaxHHProfitIncome
integer :: AdjGovBudgetConstraint,ConvergenceRelToOutput,AdjustProdGridFrisch
integer :: GovExpConstantFracOutput,GovBondResidualZeroWorld,PinKappa1ByKappa02,ImposeMaxHours,OneAssetNoCapital
integer :: ReadEarningsProcess,PerfectAnnuityMarkets,CalibrateRhoAtInitialGuess
integer :: MatchRelativeToTargetOutput,DividendFundLumpSum
logical :: stickytransition
integer :: maxiter,maxitertransflex,maxitertranssticky,maxiterKFE,maxiterequmss,maxiterrho
real(8) :: Vtol,KFEtol,deltass,deltakfe,deltacumcon,nendtrans,deltatransvec(Ttransition),cumdeltatrans(Ttransition),toltransition,deltatransparam,deltatransmin,deltatransmax,adjfricshgridfrac
real(8) :: tolequmss,stepequmss,tolrho,stepstickytransK,stepstickytransB,stepstickytransRb,dVamin,dVbmin

!CALIBRATION OPTIONS
integer ::  EstimateKappa0, EstimateKappa1, EstimateKappa2, EstimateKappa3, EstimateKappa4, EstimateRho, EstimateBorrWedge,EstimateGamma
integer :: MatchMeanIll, MatchMedianIll, MatchP75Ill, MatchFracIll0, MatchMeanLiq, MatchFracIll0Liq0, MatchMedianLiq, MatchFracLiq0, MatchFracLiqNeg,MatchKYratio
logical :: calibrating,iteratingtransition
real(8) :: defnaclose,defnbclose

!SHOCK GLOBALS
integer :: IncludeMonetaryShock
integer :: IncludeForwardGuideShock,ForwardGuideShockQtrs
real(8) :: MonetaryShockSize,ForwardGuideShockSize
real(8) :: MonetaryShockPers,ForwardGuideShockPers
logical :: forwardguide

!GRIDS GLOBALS
real(8), dimension(ngpy) 		:: ygrid,logygrid		!individual productivity
real(8), dimension(ngpa)  		:: agrid,adelta,adrift        !illiquid asset
real(8), dimension(ngpb)  		:: bgrid,bdelta,bdrift        !liquid asset
real(8), dimension(nab)  		:: abdelta
real(8), dimension(naby)  		:: abydelta
integer, dimension(naby)  		:: afromaby,bfromaby,yfromaby,abfromaby
integer, dimension(nab)  		:: afromab,bfromab
integer, dimension(ngpa,ngpb,ngpy)  :: abyfromaby
integer, dimension(ngpa,ngpb)  		:: abfromab
real(8), dimension(ngpb-1)  		:: dbgrid        !liquid asset spacing
real(8), dimension(ngpa-1)  		:: dagrid        !illiquid asset spacing

!GLOBALS FOR INCOME RISK
real(8), dimension(ngpy)		    :: ydist
real(8), dimension(ngpy,ngpy)		:: ytrans,ymarkov,ymarkovdiag,ymarkovoff

!GLOBALS FOR VALUE FUNCTIONS AND DECISION
real(8), dimension(:,:,:), allocatable		:: V,Vnew,u,gjoint,c,h,d,s,ccum1,ccum4,ccum2,dcum1,dcum4,dcum2,bdot,mpc,subeff1ass,subeff2ass,wealtheff1ass,wealtheff2ass
real(8), dimension(:,:), allocatable		:: gamarg,gbmarg,gmat,gabcum,gabmarg
real(8), dimension(:), allocatable			:: gvec

!ITERATION GLOBALS
real(8) :: delta

!PARAMETER GLOBALS
real(8)     :: rho,gam,utilcost,chi,frisch,blim,nbl,abl,prefshock,fundlev,fundbond,deathrate,meanlabeff,profdistfrac
real(8)     :: elast,alpha,deprec,alphatilde,theta,phitaylor,phifg,bondelast,borrwedge,mpshock,bondadjust,bondelastrelgdp
real(8) 	:: kappa0_d,kappa1_d,kappa2_d,kappa0_w,kappa1_w,kappa2_w,kappa3,kappa4_d,kappa4_w
real(8) 	:: dmin,taxincrstart,taxincrdecay,utilelast,utilelastalpha

!EQUILIBRIUM GLOBALS
real(8)     :: ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,profit,dividend,divrate,priceadjust,equity
real(8)     :: labtax,lumptransfer,lumptransferpc,ssdebttogdp,corptax,illassetdrop,caputil,tfpadj
integer 	:: neqmiter
logical		:: converged,initialSS

!STATISTICS GLOBALS
real(8) 	:: Ea,Eb,Ec,Elabor,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Enw,EbN,EbP,Eadjcost
real(8) 	:: FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close
real(8) 	:: PERCa(11),PERCb(11),PERCnw(11),PERCc(11),PERCinc(11)
real(8) 	:: GINIa,GINIb,GINInw,GINIc,GINIinc
real(8) 	:: Ea_nwQ(4),Eb_nwQ(4),Ec_nwQ(4),Einc_nwQ(4),Ea_incQ(4),Eb_incQ(4),Ec_incQ(4),Einc_incQ(4)
real(8)		:: Ec_bN,Ec_b0close,Ec_b0far

!CALIBRATION GLOBALS
integer					:: nparam,nmoments,objeval,ndfls
real(8), allocatable	:: paramguess(:),paramout(:),paramscale(:),paramlb(:),paramub(:),diagweight(:),nobsvec(:)
real(8)		:: targetMeanIll,targetMeanLiq,targetMedianIll,targetP75Ill,targetMedianLiq,targetFracIll0,targetFracLiq0,targetFracIll0Liq0,targetFracLiqNEG,targetKYratio
real(8)		:: modelMeanIll,modelMeanLiq,modelMedianIll,modelP75Ill,modelMedianLiq,modelFracIll0,modelFracLiq0,modelFracIll0Liq0,modelFracLiqNEG,modelKYratio
real(8) 	:: kappa2min,borrwedgemax



!SPARSE MATRIX TYPES
type COO
	integer	:: nz
	real(8), dimension(:), allocatable :: val
	integer, dimension(:), allocatable :: row, col
end type

type CSR
	integer	:: n,nz
	real(8), allocatable :: val(:)
	integer, allocatable :: row(:), col(:)
end type

!SPARSE MATRICES
type(COO), dimension(ngpy)	:: ACOO,AUCOO
type(CSR), dimension(ngpy)	:: ACSR,BCSR,AUCSR

!SOLUTION TYPE
type SolutionType
	real(8), dimension(:,:,:), allocatable	:: V,c,s,h,d,u,gjoint,bdot,mpc,subeff1ass,subeff2ass,wealtheff1ass,wealtheff2ass
	real(8), dimension(:), allocatable 		:: gvec	
	real(8), dimension(:,:), allocatable 	:: gamarg,gbmarg,gmat
	type(CSR), dimension(:), allocatable 	:: A,B,AU
end type

type(SolutionType)		:: solnINITSS,solnFINALSS,solnTRANS(Ttransition)


!EQUILIBRIUM TYPE
type EquilibriumType
	real(8)		:: ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,labtax,borrwedge,rho,mpshock,prefshock,&
					priceadjust,fundlev,elast,gam,fundbond,profit,dividend, divrate,lumptransfer,equity,caputil,deprec,tfpadj,illassetdrop
end type

type(EquilibriumType)	:: equmINITSS,equmFINALSS,equmTRANS(Ttransition)


!DISTRIBUTION STATISTICS TYPE
type DistributionStatsType
	real(8) 	:: Ea,Eb,Ec,Elabor,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, &
					EbN,EbP,Eadjcost,PERCa(11),PERCb(11),PERCnw(11),PERCc(11),PERCinc(11),GINIa,GINIb,GINInw,GINIc,GINIinc, &
					Ea_nwQ(4),Eb_nwQ(4),Ec_nwQ(4),Einc_nwQ(4),Ea_incQ(4),Eb_incQ(4),Ec_incQ(4),Einc_incQ(4),Ec_bN,Ec_b0close,Ec_b0far	
end type

type(DistributionStatsType)	::statsINITSS,statsFINALSS,statsTRANS(Ttransition)

!CUMULATIVE POLICY FUNCTION TYPE
type CumulativePolicyType
	real(8), dimension(:,:,:), allocatable		:: ccum1,ccum2,ccum4,dcum1,dcum2,dcum4
end type	

type(CumulativePolicyType)		:: cumINITSS


!IMPULSE RESPONSE FUNCTION TYPE
type ImpulseResponseType
	type(SolutionType)			:: solnSTICKY(Ttransition)
	type(EquilibriumType)		:: equmSTICKY(Ttransition)
	type(DistributionStatsType)	:: statsSTICKY(Ttransition)
	type(CumulativePolicyType)	:: cumSTICKY
end type	

type(ImpulseResponseType), dimension(0:0), target 	:: irfstruct
type(ImpulseResponseType), target					:: irfpriceexp
type(ImpulseResponseType), pointer 					:: irfsave,irfpointer,irfpointer_fs 


!THREADPRIVATE GLOBALS
real(8)							:: gbdrift,gnetwage,gill,gVb,gidioprod
!$OMP THREADPRIVATE(gbdrift,gnetwage,gill,gVb,gidioprod)

END MODULE Globals