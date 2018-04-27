%% Parameters
global LaborSupply chi frisch prefshock ScaleDisutilityIdio ImposeMaxHours ngpy ngpb ngpa nab Ttransition
global kappa0_w kappa1_w kappa2_w kappa0_d kappa1_d kappa2_d kappa3 dmax Display
global MonetaryShockPers rho rb
global dVamin dVbmin DistributeProfitsInProportion TaxHHProfitIncome profdistfrac corptax deathrate ReportNonMonotonicity
global SolveStickyPriceTransition SaveTime1PolicyFns SaveCumPolicyFnsIRF
% OPTIONS
TwoPointWageProcess = false; % set ngpy=2
Borrowing = true;

% GRID SIZES
ngpa = 40; % grid for illiquid assets
ngpbPOS = 40; % grid for liquid assets, >=0 range
ngpbNEG = 10; % grid for liquid assets, <0 range only active if Borrowing==1
ngpb = ngpbPOS + Borrowing*ngpbNEG;
ngpy = 33;
naby = ngpa*ngpb*ngpy;
nab = ngpa*ngpb;

% PARAMETERS FOR GRID CONSTRUCTION
agridparam = 0.15; % for a: approaches linear as goes to 1, approaches L shaped as goes to 0
bgridparam = 0.35; % for b pos: approaches linear as goes to 1, approaches L shaped as goes to 0
bgridparamNEG = 0.4; % for b neg: approaches linear as goes to 1, approaches L shaped as goes to 0
amax = 2000; % multiple of quarterly output
bmax = 40;

% OTHER PARAMETERS
cmin = 1e-5; % minimum consumption for natural borrowing limit
dmax = 1e10; % maximum deposit rate, for numerical stability while converging
facc = 1e-10;

Ttransition = 200; % no. time steps for the transition (each step is can be a different number of time units)

%% Globals
% OUTPUT DIR
OutputDir = '../FortranOutputDir/BaselineOutputSubdir/';
EarningsProcessDir = '../Fortran/earnings_input';
mkdir(OutputDir)

%% OPTIONS
% display options
Display = 1;
ReportNonMonotonicity = false;

% run options
CalibrateDiscountRate = false;
EquilibriumR = true;
ComputeCumulativeMPC = true;
DoImpulseResponses = true;
DoPriceExperiments = true;
SaveTime1PolicyFns = true;
SaveCumPolicyFnsIRF = false;
ComputeDiscountedMPC = true;

% labor risk options
ReadEarningsProcess = true;
LaborSupply = 1; % labor supply options: 0 = No labor supply, 1 = Separable, 2 = GHH
ScaleDisutilityIdio = false;
ImposeMaxHours = true;
PerfectAnnuityMarkets = true;
GovBondResidualZeroWorld= true; % imposes closed economy and imputes residual bond holdings to govt
AdjustProdGridFrisch = true;
adjfricshgridfrac = 0.85; % fraction of Frisch to adjust by

% calibration options
CalibrateCostFunction = false;
CalibrateRhoAtInitialGuess = false;
MatchRelativeToTargetOutput = false;
ImposeEqumInCalibration = false;

% adjustment cost options
OneAssetNoCapital = false;
PinKappa1ByKappa02 = true;

% transition computation options
SolveStickyPriceTransition = true;
ConvergenceRelToOutput = true;
FirmDiscountRate = 5; % 1 for rho, 2 for rb initial steady state, 3 for ra initial steady state, 4 for rb transition, 5 for ra transition
bondelastrelgdp = 1; % bigger for smaller interest rate movements, closer to zero for larger interest rate movements. relative to steady state gdp
bondadjust = 0.1; % more responsive interest rate when closer to zero

% dividend options
DividendFundLumpSum = true;
DistributeProfitsInProportion = true; % distributes profdistfrac of profits as dividends, and (1-profdistfrac) to households in proportion to productivity
profdistfrac = 0.33; % set to alpha to neutralize effect of profits on redistribution between liquid and illiquid assets
TaxHHProfitIncome = true; % taxes profit income at labor tax rate if DistributeProfitsInProportion = true


% government bc options
AdjGovBudgetConstraint = 2; % 1 for adjust spending, 2 for adjust lump sum taxes, 3 for let debt adjust (choose options below for financing), 4 for adjust proportional tax
GovExpConstantFracOutput = false; % only active if AdjGovBudgetConstraint==3
taxincrstart = 1; % quarters after shock that fiscal policy adjusts
taxincrdecay = 0.02; % decay rate for tax increase higher for faster decay


% MONETARY POLICY SHOCK
IncludeMonetaryShock = true;
MonetaryShockSize = -0.0025; % percentage points
MonetaryShockPers = exp(-0.5); % quarterly

% FORWARD GUIDANCE SHOCK
IncludeForwardGuideShock = false;
ForwardGuideShockSize = -0.0025; % percentage points
ForwardGuideShockPers = exp(-0.5); % quarterly
ForwardGuideShockQtrs = 9; % number of quarters in advance (set phifg below)


% CALIBRATION OPTIONS
EstimateKappa0 = true;
EstimateKappa1 = false;
EstimateKappa2 = true;
EstimateKappa3 = false;
EstimateKappa4 = false;
EstimateRho = true;
EstimateBorrWedge = true;
EstimateGamma = false;


MatchMeanIll = false;
MatchKYratio = true;
MatchMedianIll = false;
MatchP75Ill = false;
MatchFracIll0 = false;
MatchMeanLiq = true;
MatchMedianLiq = false;
MatchFracLiq0 = true;
MatchFracLiqNeg = true;
MatchFracIll0Liq0 = true;

defnbclose = 0;
defnaclose = 0;

ndfls = 2;

% SOLUTION PARAMETERS
maxiter = 500;
Vtol = 1e-8;
maxiterKFE = 2000;
KFEtol = 1e-12;
deltass = 1e6;
deltakfe = 1e6;
dVamin = 1e-8;
dVbmin = 1e-8;

tolequmss = 1e-10;
stepequmss = 0.1;
maxiterequmss = 40;

toltransition = 1e-6;
deltatransmin = 1/3;
deltatransmax = 40;
deltatransparam = 0.35;
maxitertranssticky = 5000;
stepstickytransK = 0.01;
stepstickytransB = 0.001;

deltacumcon = 0.01; % set to a low number like 0.01 for accurate steady state MPCs, and to deltratransmin for IRF consistency

% discount rates
rho = 0.01272892513; % (baseline), 0.0133348071(frisch=0.5), 0.0123887061(frisch=1.5), 0.0125456234 (ghh)

% preferences
deathrate = 1/(4*45); % poisson death rate
gam = 1; % risk aversion
prefshock = 1;

% liquid assets
rb = 0.02/4; % liquid return
borrwedge = 0.0148846; % quarterly wedge between rb and rborr: intermediation cost
borrwedgemax = 0.09;
blim = -1; % borrowing limit multiple of quarterly output

rborr = rb + borrwedge;

% withdrawal costs
kappa0_w = 0.04383;
kappa2_w = 0.40176;
kappa3 = 0.03*2.92/4; % as a fraction of average illiquid assets
kappa4_w = 0;

if PinKappa1ByKappa02
	kappa1_w = ((1-kappa0_w)*(1+kappa2_w))^(-1/kappa2_w);
else
    kappa1_w = 0.48236; %#ok<*UNRCH>
end

kappa2min = 0.05; % to make sure there is enough curvature for calibration

% deposit costs
kappa0_d = kappa0_w;
kappa1_d = kappa1_w;
kappa2_d = kappa2_w;
kappa4_d = kappa4_w;

theta = 100; % price adjustment cost
phitaylor = 1.25; % inflation coefficient in taylor rule, 0 for fixed nominal interest rate at ss.
phifg = phitaylor; % with forward guidance use phifg instead of phitaylor pre-shock, and phitaylor after (not equal 1)

pi = 0; % inflation rate: this is steady state target 
rnom = rb + pi; % nominal interest rate (fisher equation): this will be constant in taylor rule
mpshock = 0; % shock to nominal interest rate in taylor rule

tfp = 1;
elast = 10; % elasticity of DS aggregator
gap = 0; % steady state output gap
mc = (elast-1)/elast;

alpha = 0.33;
alphatilde = (alpha^alpha)*((1-alpha)^(1-alpha));
deprec = 0.07/4; % depreciation rate

fundlev = 0;
utilelast = 0;
utilelastalpha = 1 + utilelast-alpha*utilelast;

% government
labtax = 0.30;
lumptransferpc = 40000*labtax/(115000*2.92*4);
corptax = 0;
ssdebttogdp = 0.26*4; % if foreign sector assumed to hold residual bonds

% calibration targets (relative to quarterly output)
targetMeanIll = 2.92 * 4;
targetMeanLiq = 0.2 * 4;
targetMedianIll = 0.21 * targetMeanIll;
targetP75Ill = 0.71 * targetMeanIll;
targetMedianLiq = 0.085 * targetMeanLiq;
targetFracIll0 = 0.115 + 0.12;
targetFracLiq0 = 0.35;
targetFracLiqNEG = 0.15;
targetFracIll0Liq0 = 0.10;

lumptransfer = lumptransferpc*targetMeanIll;

if DividendFundLumpSum
    la = -(deprec+rb*fundlev);
    lb = (elast-1)/elast*alpha + targetMeanIll/(1-fundlev)*(deprec+rb*fundlev);
    if DistributeProfitsInProportion
        lb = lb + (1/elast)*(1-corptax)*profdistfrac;
    else
        lb = lb + (1/elast)*(1-corptax);
    end
    lc = -targetMeanIll/(1-fundlev)*(elast-1)/elast*alpha;
    targetKYratio = (-lb+sqrt(lb^2-4*la*lc))/(2*la);
else
    targetKYratio = targetMeanIll/(1-fundlev);
end

% if solving for equilibrium, these are guesses
KYratio = targetKYratio;
KNratio = (tfp*KYratio)^(1/(1-alpha));
rcapital = mc*alpha/KYratio;
wage = mc*(1-alpha)*tfp*(KNratio^alpha);
netwage = (1-labtax)*wage;
if DividendFundLumpSum
    divrate = 0;
else
    divrate = (1-corptax)*(1-mc)/KYratio; % outside of steady state include price adjustments
end
if DistributeProfitsInProportion
    divrate = profdistfrac*divrate;
end

ra = (rcapital - deprec + divrate - fundlev*rb)/(1-fundlev);

if LaborSupply
    meanlabeff = (KYratio/KNratio)*3; % scale efficiency units so that average hours would be 1/3 in eqm if cov(z,h)=0
else
    meanlabeff = KYratio^(alpha/(alpha-1)); % scale efficiency units so that output equals 1
end

frisch = 1; % frisch elasticity labor supply

if OneAssetNoCapital
    alpha = 0;
    profdistfrac = 0;
    rcapital = 0;
    wage = mc*(1-alpha)*tfp;
    netwage = (1-labtax)*wage;
    ra = 0;
    meanlabeff = 3;
    kappa0_w = 100;
    kappa0_d = 100;
    FirmDiscountRate = 4;
end

% guess chi's so that at average wages and average consumption hours =1/3 (sets C/Y = 0.75)
switch LaborSupply
    case 0; chi = 0; % No
    case 1; chi = meanlabeff * 0.75^gam * 3^(1/frisch); % Sep
    case 2; chi = meanlabeff * 3^(1/frisch); % GHH
end

% requires an intial guess of profits if DistributeProfitsInProportion = 1, since output not known in advance
if DistributeProfitsInProportion
    if OneAssetNoCapital
        profit = 1-mc;
    else
        % these will be updated in iterations
        switch LaborSupply
            case 0; capital = 12; % No
            case 1; capital = 16; % Sep
            case 2; capital = 14; % GHH
        end
        profit = (1-mc)*capital/KYratio;
    end
end

calibrating = false;
iteratingtransition = false;
