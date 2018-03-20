if Display; disp('Solving for initial steady state'); end
initialSS = true;

if CalibrateDiscountRate
    neqmiter = 1;
    fopen(fullfile(OutputDir,'DiscountRateCalibration.txt'),'w');
    if ngpy==1 && deathrate==0.0
        % RA model no death
        lrho = icdf('Logistic',exp([-0.02 -0.015]));
        rho = fzero(@FnDiscountRate,lrho,optimset('TolX',1e-6));
    else
        lrho = icdf('Logistic',exp([-0.016 -0.014]));
        rho = fzero(@FnDiscountRate,lrho,optimset('TolX',1e-8));
    end
end
if EquilibriumR
    SolveSteadyStateEqum
else
    Grids
    IterateBellman
    StationaryDistribution
end

if ComputeCumulativeMPC; CumulativeConsumption; end
if ComputeDiscountedMPC
    solnINITSS = DiscountedMPC(c,d,AU);
end

% solnINITSS.A = A;
solnINITSS.AU = AU;
solnINITSS.V = V;
% solnINITSS.u = u;
solnINITSS.c = c;
solnINITSS.s = s;
solnINITSS.h = h;
solnINITSS.d = d;
solnINITSS.bdot = bdot;

% solnINITSS.gjoint = gjoint;
solnINITSS.gamarg = gamarg;
solnINITSS.gbmarg = gbmarg;
% solnINITSS.gvec = gvec;
solnINITSS.gmat = gmat;
solnINITSS.gabmarg = gabmarg;
solnINITSS.gabcum = gabcum;

DistributionStatistics
% EquilibriumType
equmINITSS = v2struct( ra,rborr,rcapital,wage,netwage,KYratio,KNratio,mc,rb,tfp,pi,rnom,gap,bond,capital,labor,output,investment,govexp,taxrev,govbond,worldbond,labtax, ...
        borrwedge,rho,mpshock,prefshock,priceadjust,fundlev,elast,gam,fundbond,profit,dividend,divrate,lumptransfer,equity,caputil,deprec,tfpadj,illassetdrop);
% DistributionStatsType
statsINITSS = v2struct(Ea,Eb,Ec,Elabor,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, ...
         EbN,EbP,Eadjcost,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, ...
         Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far);

initialSS = false;
