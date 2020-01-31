if iteratingtransition
    equmTRANS.rborr = equmTRANS.rb + equmTRANS.borrwedge;

    % world bond
    equmTRANS.worldbond = [equmINITSS.worldbond WorldBondFunction2(equmTRANS.rb(1:Ttransition-1),equmINITSS.worldbond,equmINITSS.rb,bondelast)];
    for it = 1:Ttransition-1
        equmTRANS.worldbond(it+1) = equmTRANS.worldbond(it) + bondadjust*deltatransvec(it)*(equmTRANS.worldbond(it+1)-equmTRANS.worldbond(it));
    end

    % fund bond
    equmTRANS.fundbond = -equmTRANS.capital.*equmTRANS.fundlev;

    % solve phillips curve backwards for marginal costs
    switch FirmDiscountRate
        case 1; lfirmdiscount = equmTRANS.rho;
        case 2; lfirmdiscount = equmINITSS.rb;
        case 3; lfirmdiscount = equmINITSS.ra;
        case 4; lfirmdiscount = equmTRANS.rb;
        case 5; lfirmdiscount = equmTRANS.ra;
    end

    % marginal costs
    equmTRANS.mc = max(lminmargcost,(lfirmdiscount-diff([equmTRANS.tfp equmFINALSS.tfp])./(equmTRANS.tfp.*deltatransvec) ...
        - alpha*diff([equmTRANS.capital equmFINALSS.capital])./(equmTRANS.capital.*deltatransvec) ...
        - alpha*diff([equmTRANS.caputil equmFINALSS.caputil])./(equmTRANS.caputil.*deltatransvec) ...
        - (1-alpha)*diff([equmTRANS.labor equmFINALSS.labor])./(equmTRANS.labor.*deltatransvec) ) ...
        .* [equmTRANS.pi(2:Ttransition) equmFINALSS.pi] * theta./ equmTRANS.elast ...
        + (equmTRANS.elast-1)./equmTRANS.elast - (diff([equmTRANS.pi equmFINALSS.pi])./deltatransvec) * theta./ equmTRANS.elast);

    equmTRANS.gap = equmTRANS.elast.*equmTRANS.mc ./ (equmTRANS.elast-1) - 1;
    equmTRANS.tfpadj = (equmTRANS.tfp.^((1+utilelast)/utilelastalpha)) .* (equmTRANS.mc*alpha./equmINITSS.rcapital).^(alpha*utilelast/utilelastalpha);
    equmTRANS.KNratio = equmTRANS.capital./equmTRANS.labor;
    equmTRANS.wage = equmTRANS.mc*(1-alpha).* equmTRANS.tfpadj .* equmTRANS.KNratio.^(alpha/utilelastalpha);
    equmTRANS.netwage = (1-equmTRANS.labtax).*equmTRANS.wage;
    equmTRANS.caputil = ((equmTRANS.mc*alpha.*equmTRANS.tfp/equmINITSS.rcapital) .* equmTRANS.KNratio.^(alpha-1)) .^ (utilelast/utilelastalpha);
    equmTRANS.output = equmTRANS.tfpadj .* equmTRANS.capital.^(alpha/utilelastalpha) .* equmTRANS.labor.^((1-alpha)*(1+utilelast)/utilelastalpha);
    equmTRANS.KYratio = (equmTRANS.KNratio.^(1-alpha)) ./ (equmTRANS.tfp.*equmTRANS.caputil.^alpha);
    equmTRANS.rcapital = (equmINITSS.rcapital.^utilelast * equmTRANS.mc * alpha ./ equmTRANS.KYratio) .^ (1/(1+utilelast));
    equmTRANS.priceadjust = theta/2*equmTRANS.pi.^2.*equmTRANS.capital./equmTRANS.KYratio;
    equmTRANS.profit = (1-equmTRANS.mc).*equmTRANS.capital./equmTRANS.KYratio - equmTRANS.priceadjust;

    equmTRANS.deprec = equmINITSS.deprec + (utilelast*equmINITSS.rcapital/(1+ utilelast)) .* ((equmTRANS.rcapital/equmINITSS.rcapital).^(1+utilelast) -1);

    % solve backward for investment
    equmTRANS.investment = diff([equmTRANS.capital equmFINALSS.capital])./deltatransvec + equmTRANS.deprec.*equmTRANS.capital;

    % dividends and illiquid return
    equmTRANS.dividend = equmTRANS.profit*(1-corptax);
    if DistributeProfitsInProportion; equmTRANS.dividend = profdistfrac*equmTRANS.dividend; end

    if DividendFundLumpSum; equmTRANS.divrate = 0;
    else; equmTRANS.divrate = equmTRANS.dividend./equmTRANS.capital; end

    equmTRANS.ra = (equmTRANS.rcapital.*equmTRANS.caputil - equmTRANS.deprec + equmTRANS.divrate - equmTRANS.fundlev.*equmTRANS.rb) ./ (1-equmTRANS.fundlev);

    % value of equity component of investmemt fund
    if DividendFundLumpSum
        it = Ttransition;
        equmTRANS.equity(it) = (equmFINALSS.equity + equmTRANS.dividend(it)*deltatransvec(it)) / (1+deltatransvec(it)*equmTRANS.ra(it));
        for it = Ttransition-1:-1:1
            equmTRANS.equity(it) = (equmTRANS.equity(it+1) + equmTRANS.dividend(it)*deltatransvec(it)) / (1+deltatransvec(it)*equmTRANS.ra(it));
        end
        equmTRANS.illassetdrop = ((1-equmTRANS.fundlev(1))*equmTRANS.capital(1) + equmTRANS.equity(1)) / ((1-equmINITSS.fundlev)*equmINITSS.capital + equmINITSS.equity);
    else
        equmTRANS.equity = 0;
        equmTRANS.illassetdrop = 1;
    end

    % government budget constraint,expenditures and tax rates
    switch AdjGovBudgetConstraint
        case 1 % adjust spending
            equmTRANS.govbond = equmINITSS.govbond;
            equmTRANS.labtax = equmINITSS.labtax;
            equmTRANS.lumptransfer = equmINITSS.lumptransfer;

            equmTRANS.taxrev = equmTRANS.labtax*equmTRANS.wage*equmTRANS.labor - equmTRANS.lumptransfer + corptax*equmTRANS.profit;
            if DistributeProfitsInProportion && TaxHHProfitIncome; equmTRANS.taxrev = equmTRANS.taxrev + equmTRANS.labtax*(1-profdistfrac)*equmTRANS.profit*(1-corptax); end

            equmTRANS.govexp = equmTRANS.taxrev + equmTRANS.rb*equmINITSS.govbond;

        case 2 % adjust lump sum taxes
            equmTRANS.govbond = equmINITSS.govbond;
            equmTRANS.govexp = equmINITSS.govexp;
            equmTRANS.labtax = equmINITSS.labtax;
            equmTRANS.taxrev = equmTRANS.govexp - equmTRANS.rb*equmINITSS.govbond;
            equmTRANS.lumptransfer = equmTRANS.labtax*equmTRANS.wage.*equmTRANS.labor + corptax*equmTRANS.profit + equmTRANS.rb*equmINITSS.govbond - equmTRANS.govexp;
            if DistributeProfitsInProportion && TaxHHProfitIncome; equmTRANS.lumptransfer = equmTRANS.lumptransfer + equmTRANS.labtax*(1-profdistfrac)*equmTRANS.profit*(1-corptax); end

        case 3 % adjust debt
            if GovExpConstantFracOutput; equmTRANS.govexp = equmTRANS.output*equmINITSS.govexp/equmINITSS.output;
            else; equmTRANS.govexp = repmat(equmINITSS.govexp,1,Ttransition); end

            equmTRANS.lumptransfer = equmINITSS.lumptransfer;
            equmTRANS.taxrev = equmTRANS.labtax*equmTRANS.wage*equmTRANS.labor - equmTRANS.lumptransfer + corptax*equmTRANS.profit;
            if DistributeProfitsInProportion && TaxHHProfitIncome; equmTRANS.taxrev = equmTRANS.taxrev + equmTRANS.labtax*(1-profdistfrac)*equmTRANS.profit*(1-corptax); end

            % compute required increase in lumptransfer
            lrgov = equmTRANS.rb;
            lpvgovbc = equmFINALSS.govbond;
            lpvlumpincr = 0;
            for it = Ttransition:-1:1
                lpvgovbc = (lpvgovbc + deltatransvec(it)*(equmTRANS.govexp(it) - equmTRANS.taxrev(it)))/(1+deltatransvec(it)*lrgov(it));
                if cumdeltatrans(it)<taxincrstart; lpvlumpincr = lpvlumpincr/(1+deltatransvec(it)*lrgov(it));
                else; lpvlumpincr = (lpvlumpincr + deltatransvec(it))/(1+deltatransvec(it)*(lrgov(it)+taxincrdecay)); end
            end

            linitlumpincr = (equmINITSS.govbond-lpvgovbc) / lpvlumpincr;
            equmTRANS.lumptransfer = equmINITSS.lumptransfer + linitlumpincr*exp(-taxincrdecay*(cumdeltatrans-taxincrstart));
            equmTRANS.lumptransfer(cumdeltatrans<taxincrstart) = equmINITSS.lumptransfer;

            equmTRANS.taxrev = equmTRANS.labtax*equmTRANS.wage*equmTRANS.labor - equmTRANS.lumptransfer + corptax*equmTRANS.profit;
            if DistributeProfitsInProportion && TaxHHProfitIncome; equmTRANS.taxrev = equmTRANS.taxrev + equmTRANS.labtax*(1-profdistfrac)*equmTRANS.profit*(1-corptax); end

            equmTRANS.govbond(Ttransition) = equmFINALSS.govbond;
            for it = Ttransition-1:-1:2
                equmTRANS.govbond(it) = (equmTRANS.govbond(it+1) - deltatransvec(it)*(equmTRANS.taxrev(it)-equmTRANS.govexp(it))) / (1+deltatransvec(it)*lrgov(it));
            end
            equmTRANS.govbond(1) = equmINITSS.govbond;

            equmTRANS.lumptransfer = equmTRANS.lumptransfer + (equmTRANS.rb-lrgov).*equmTRANS.govbond;
            equmTRANS.taxrev = equmTRANS.labtax*equmTRANS.wage*equmTRANS.labor - equmTRANS.lumptransfer + corptax*equmTRANS.profit;
            if DistributeProfitsInProportion && TaxHHProfitIncome; equmTRANS.taxrev = equmTRANS.taxrev + equmTRANS.labtax*(1-profdistfrac)*equmTRANS.profit*(1-corptax); end

        case 4 % adjust proportional tax rate
            equmTRANS.govbond = equmINITSS.govbond;
            equmTRANS.govexp = equmINITSS.govexp;
            equmTRANS.lumptransfer = equmINITSS.lumptransfer;
            equmTRANS.taxrev = equmTRANS.govexp - equmTRANS.rb*equmINITSS.govbond;

            if DistributeProfitsInProportion && TaxHHProfitIncome; equmTRANS.labtax  = (equmTRANS.lumptransfer - corptax*equmTRANS.profit - equmTRANS.rb*equmINITSS.govbond + equmTRANS.govexp) / (equmTRANS.wage*equmTRANS.labor + (1-profdistfrac)*equmTRANS.profit*(1-corptax));
            else; equmTRANS.labtax  = (equmTRANS.lumptransfer - corptax*equmTRANS.profit - equmTRANS.rb*equmINITSS.govbond + equmTRANS.govexp) / (equmTRANS.wage*equmTRANS.labor); end
    end

    % household bonds
    equmTRANS.bond = -equmTRANS.worldbond - equmTRANS.govbond - equmTRANS.fundbond;

    holdbgrid = bgrid;

    %% check Fortran
    TRANS = load(sprintf('TRANS/sol%4d.txt',ii));
    TRANS = struct('borrwedge',TRANS(1,:),'fundlev',TRANS(2,:),'elast',TRANS(3,:),'tfp',TRANS(4,:),'caputil',TRANS(5,:),'labor',TRANS(6,:),'labtax',TRANS(7,:),'ra',TRANS(8,:),'mpshock',TRANS(9,:),'capital',TRANS(10,:),'pi',TRANS(11,:),'rnom',TRANS(12,:),'rb',TRANS(13,:),'rborr',TRANS(14,:),'worldbond',TRANS(15,:),'fundbond',TRANS(16,:),'mc',TRANS(17,:),'gap',TRANS(18,:),'tfpadj',TRANS(19,:),'KNratio',TRANS(20,:),'wage',TRANS(21,:),'netwage',TRANS(22,:),'output',TRANS(23,:),'KYratio',TRANS(24,:),'rcapital',TRANS(25,:),'priceadjust',TRANS(26,:),'profit',TRANS(27,:),'deprec',TRANS(28,:),'investment',TRANS(29,:),'dividend',TRANS(30,:),'divrate',TRANS(31,:),'equity',TRANS(32,:),'illassetdrop',TRANS(33,:),'govbond',TRANS(34,:),'govexp',TRANS(35,:),'taxrev',TRANS(36,:),'lumptransfer',TRANS(37,:),'bond',TRANS(38,:));
    for f=fields(TRANS)'
        err.(f{1}) = max(abs(TRANS.(f{1})-equmTRANS.(f{1})));
    end
    disp(err)
%     plot(equmTRANS.worldbond-TRANS.worldbond)
end

%% solve backward
B = nan(nab,5,ngpy,Ttransition);
for it = Ttransition:-1:1
    if Display>1; fprintf('  Solving transition backward: %d\n',it); end
    if it==Ttransition; V = solnFINALSS.V; else; V = solnTRANS.V{it+1}; end

    % set drifts and globals
%     rho = equmTRANS.rho(it);
    ra = equmTRANS.ra(it);
    rborr = equmTRANS.rborr(it);
    borrwedge = equmTRANS.borrwedge;
    wage = equmTRANS.wage(it);
    netwage = equmTRANS.netwage(it);
    labtax = equmTRANS.labtax;
    lumptransfer = equmTRANS.lumptransfer(it);
    rb = equmTRANS.rb(it);
    tfp = equmTRANS.tfp;
    mpshock = equmTRANS.mpshock(it);
%     prefshock = equmTRANS.prefshock(it);
    fundlev = equmTRANS.fundlev;
%     fundbond = equmTRANS.fundbond;
    worldbond = equmTRANS.worldbond(it);
    elast = equmTRANS.elast;
%     gam = equmTRANS.gam(it);
    profit = equmTRANS.profit(it);

    % set drifts
    ltau = 15;
    ltau0 = (ra+PerfectAnnuityMarkets*deathrate)*(amax*0.999)^(1-ltau);
    adrift = (ra+PerfectAnnuityMarkets*deathrate)*agrid - ltau0*agrid.^ltau;
    bdrift = (rborr*(bgrid<=0)+rb*(bgrid>0)+PerfectAnnuityMarkets*deathrate).*bgrid;

    nblviolated = false;
    lbgrid{it} = bgrid;
    if Borrowing && rborr+PerfectAnnuityMarkets*deathrate>0 && bgrid(1) < -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)
        if Display
            fprintf('Warning: natural borrowing limit violated in transition\n')
            fprintf('Steady State ABL: %.15g   Current NBL: %15.g\n',abl,-lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate))
            fprintf(' pi %.15g\n',equmTRANS.pi(it))
            fprintf(' rborr %.15g\n',rborr)
            fprintf(' rb %.15g\n',rb)
            fprintf(' wage %.15g\n',wage)
            fprintf(' mc %.15g\n',equmTRANS.mc(it))
            fprintf(' KYratio %.15g\n',equmTRANS.KYratio(it))
            fprintf(' mpshock %.15g\n',equmTRANS.mpshock(it))
        end
        error('STOP')
        nblviolated = true;
        lbmin = -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate) + cmin;

        lbgrid{it}(1:ngpbNEG/2+1) = PowerSpacedGrid(ngpbNEG/2+1,bgridparamNEG,lbmin,lbmin/2);
        lbgrid{it}(ngpbNEG/2+2:ngpbNEG) = lbmin-fliplr(lbgrid{it}(2:ngpbNEG/2,it));
        
        bdrift = (rborr*(lbgrid{it}<=0)+rb*(lbgrid{it}>0)+PerfectAnnuityMarkets*deathrate).*lbgrid{it};
    end

    delta = deltatransvec(it);

%     max(max(max(abs(V-reshape(load(sprintf('V %3d%3d.txt',ii,it)),ngpa,ngpb,ngpy)))))
%     max(abs(adrift-load(sprintf('adrift %3d%3d.txt',ii,it))))
%     max(abs(bdrift'-load(sprintf('bdrift %3d%3d.txt',ii,it))))
%     max(abs([netwage,lumptransfer,labtax,profit,meanlabeff,delta,ra]-load(sprintf('netwage %3d%3d.txt',ii,it))))
    [Vnew,c,h,s,d,dadj] = HJBUpdate(V,adrift,bdrift,netwage,lumptransfer,labtax,profit,meanlabeff,delta);
%     max(max(max(abs(Vnew-reshape(load(sprintf('Vnew %3d%3d.txt',ii,it)),ngpa,ngpb,ngpy)))))

    % a drift, upwind
    laudriftB = min(d + adrift,0);
    laudriftB(:,ngpb,:) = laudriftB(:,ngpb-1,:);
    laudriftF = max(d + adrift,0);
    laudriftF(:,ngpb,:) = laudriftF(:,ngpb-1,:);

    % b drift,upwind
    bdot = s - dadj;
    lbudriftB = min(bdot,0);
    lbudriftB(:,ngpb,:) = min(s(:,ngpb,:) - dadj(:,ngpb-1,:),0);
    lbudriftF = max(bdot,0);
    % lbudriftF(:,ngpb,:) = max(s(:,ngpb,:) - dadj(:,ngpb-1,:),0);

    AU = reshape([ laudriftF./[dagrid;1] ... % a+1
                  -laudriftB./[1;dagrid] ... % a-1
                   lbudriftF./[dbgrid 1] ... % b+1
                  -lbudriftB./[1 dbgrid] ... % b-1
                ], nab,4,ngpy);
    B(:,:,:,it) = [1+deltatransvec(it)*(sum(AU,2)+deathrate-ymarkovdiag) ... % diagonal value
                    -deltatransvec(it)*AU.*ABdelta]; % adjust for non-linearly spaced grids

    % store value functions and A matrix
    solnTRANS.V{it} = Vnew;
    solnTRANS.AU{it} = AU;
    solnTRANS.c{it} = c;
    solnTRANS.s{it} = s;
    solnTRANS.h{it} = h;
    solnTRANS.d{it} = d;
    solnTRANS.bdot{it} = bdot;
end

% B=[B(:,1:2,:,:) [zeros(1,1,ngpy,Ttransition);B(1:end-1,3,:,:)] B(:,4,:,:) [zeros(ngpa,1,ngpy,Ttransition);B(1:end-ngpa,5,:,:)]]; % bug?

%% simulate forward
for it = 1:Ttransition
    if Display>1; fprintf('Iterating transition forward: %d\n',it); end

    if it>1
        lgmat = solnTRANS.gmat{it-1};
    elseif ~DividendFundLumpSum || OneAssetNoCapital
        lgmat = solnINITSS.gmat;
    else
        lgmat = AdjustDistProportionately(agrid,adelta,solnINITSS.gmat,equmTRANS.illassetdrop(1));
    end

    lmat = eye(ngpy) + deltatransvec(it)*ymarkovoff';

    lgmat1 = reshape(lgmat,nab,ngpy)*lmat';
    ib = 1+ngpa*ngpbNEG;
    lgmat1(ib,:) = lgmat1(ib,:) + deltatransvec(it)*deathrate./abdelta(ib)*(abdelta(:)'*reshape(gmat,nab,ngpy));

% ACOO_ = spconvert(load('TRANS/ACOO.txt'));
% A = [-sum(A,2) A];
% ACOO = spdiags(A(:,:,1),[0 -1 1 -ngpa ngpa],nab,nab)';
% 
iba = reshape(1:nab,ngpb,ngpa)';
% max(max(abs(ACOO-ACOO_(iba,iba))))
% full(ACOO(1:5,1:5))
% full(ACOO_(iba(1:5),iba(1:5)))
% 
BCSR_ = spconvert(load('TRANS/BCSR.txt'));
BCSR = spdiags(B(:,:,1,1),[0 -1 1 -ngpa ngpa],nab,nab);
% 
max(max(abs(BCSR-BCSR_(iba,iba))))
full(BCSR(1:5,1:5))
full(BCSR_(iba(1:5),iba(1:5)))
% 
    
    % sweep over y
    for iy = 1:ngpy %par
        lgmat1(:,iy) = spdiags(B(:,:,iy,it),[0 -1 1 -ngpa ngpa],nab,nab)\lgmat1(:,iy);
    end
    lgmat1 = reshape(lgmat1,ngpa,ngpb,ngpy);
    solnTRANS.gmat{it} = lgmat1;

%     lgmat1_=permute(reshape(load(sprintf('gmat %3d%3d.txt',ii,it)),ngpb,ngpa,ngpy),[2 1 3]);
%     if any(abs(lgmat1_(39,10,:)-lgmat1(39,10,:))>0.01); lgmat1_(39,10,:)=lgmat1(39,10,:); end
%     for i=38:50
%         if any(abs(lgmat1_(40,i,:)-lgmat1(40,i,:))>0.01); lgmat1_(40,i,:)=lgmat1(40,i,:); disp(i); end
%     end
%     max(max(max(abs(lgmat1-lgmat1_))))

    % marginal distributions
    solnTRANS.gamarg{it} = sum(lgmat1.*bdelta,2);
    solnTRANS.gbmarg{it} = sum(lgmat1.*adelta,1);
    
    % set globals: apply current period policy function to last period distribution of state variables
        
    if it>1
        gmat = solnTRANS.gmat{it-1};
        gamarg = solnTRANS.gamarg{it-1};
        gbmarg = solnTRANS.gbmarg{it-1};
    elseif ~DividendFundLumpSum || OneAssetNoCapital
        gmat = solnINITSS.gmat;
        gamarg = solnINITSS.gamarg;
        gbmarg = solnINITSS.gbmarg;
    else
        gmat = lgmat;
        gamarg = sum(gmat.*bdelta,2);
        gbmarg = sum(gmat.*adelta,1);
    end    
        
    wage = equmTRANS.wage(it);
    netwage = equmTRANS.netwage(it);
    lumptransfer = equmTRANS.lumptransfer(it);
    rb = equmTRANS.rb(it);
    ra = equmTRANS.ra(it);
    rborr = equmTRANS.rborr(it);
    profit = equmTRANS.profit(it);
    
    c = solnTRANS.c{it};
    h = solnTRANS.h{it};
    d = solnTRANS.d{it};

    if it>1; bgrid = lbgrid{it-1}; end
    % PROBLEM HERE SINCE bdelta IS USED IN DISTRIBUTION STATISTICS
% if ~iteratingtransition
%     fprintf('%d\t%g\t%g\n',it,...
%         max(max(max(abs(h-reshape(load(sprintf('TRANS/h%4d.txt',it)),ngpa,ngpb,ngpy))))),...
%         max(max(max(abs(gmat-permute(reshape(load(sprintf('TRANS/gmat%4d.txt',it)),ngpb,ngpa,ngpy),[2 1 3]))))))
% end

    DistributionStatistics

    bgrid = holdbgrid;
    statsTRANS(it) = v2struct(Ea,Eb,Ec,Elabor,Ed,Ewage,Enetlabinc,Egrosslabinc,Enetprofinc,Egrossprofinc,Einc,Ehours,Enw,FRACa0,FRACa0close,FRACb0,FRACb0close,FRACb0a0,FRACb0aP,FRACbN,FRACnw0,FRACnw0close,FRACb0a0close, ...
        EbN,EbP,Eadjcost,PERCa,PERCb,PERCnw,PERCc,PERCinc,GINIa,GINIb,GINInw,GINIc,GINIinc, ...
        Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ,Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ,Ec_bN,Ec_b0close,Ec_b0far);
end