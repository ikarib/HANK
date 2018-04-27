% vectorize everything
lnw = agrid+bgrid;
ladjcost = adjcostfn(d,agrid);
lwage = ygrid*wage;
llabor = h.*ygrid;
lgrosslabinc = llabor*wage;
lnetlabinc = llabor*netwage + lumptransfer;

if DistributeProfitsInProportion
    lgrossprofinc = (1-profdistfrac)*profit*ygrid/meanlabeff;
else
    lgrossprofinc = 0;
end
if TaxHHProfitIncome
    lnetprofinc = (1-labtax)*lgrossprofinc;
else
    lnetprofinc = lgrossprofinc;
end

lgrossinc = lgrosslabinc + lgrossprofinc + bdrift + ...
    (ra+PerfectAnnuityMarkets*deathrate)*agrid;

% marginal distributions
gamarg = sum(gmat.*bdelta,2);
gbmarg = sum(gmat.*adelta,1);
lmargdist = gmat.*abdelta;

% ab joint distribution
gabmarg = sum(gmat,3).*abdelta;

% ab cumulative distribution
gabcum = cumsum(cumsum(gabmarg),2);

% liquid wealth marginal dist
lbmargdist = sum(gbmarg,3).*bdelta;
lbmargcum = cumsum(lbmargdist);

% illiquid wealth marginal dist
lamargdist = sum(gamarg,3).*adelta;
lamargcum = cumsum(lamargdist);

% other marginal dist
if ~iteratingtransition && ~calibrating
    % networth grid and marginal dist
    [lnwgrid,ordernw] = sort(lnw(:));
    lnwdelta = diff2(lnwgrid);
    lnwmargdist = gabmarg(ordernw);
    lnwmargcum = cumsum(lnwmargdist);
    [ia,ib] = ind2sub([ngpa,ngpb],ordernw);
    iab = lnwmargdist <= 1e-12;
    lnw_a = agrid(ia); lnw_a(iab) = 0;
    lnw_b = bgrid(ib)'; lnw_b(iab) = 0;
    lnw_c = sum(c.*lmargdist,3)./gabmarg;
    lnw_c = lnw_c(ordernw); lnw_c(iab) = 0;
    lnw_h = sum(h.*lmargdist,3)./gabmarg;
    lnw_h = lnw_h(ordernw); lnw_h(iab) = 0;
    lnw_inc = sum(lgrossinc.*lmargdist,3)./gabmarg;
    lnw_inc = lnw_inc(ordernw); lnw_inc(iab) = 0;

    % consumption
    [lcgrid,orderc] = sort(c(:));
    lcdelta = diff2(lcgrid);
    lcmargdist = lmargdist(orderc);
    lcmargcum = cumsum(lcmargdist);

    % gross total income
    [lincgrid,orderinc] = sort(lgrossinc(:));
    lincdelta = diff2(lincgrid);
    lincmargdist = lmargdist(orderinc);
    lincmargcum = cumsum(lincmargdist);
    [ia,ib,~] = ind2sub([ngpa,ngpb,ngpy],orderinc);
    [iab,~] = ind2sub([nab,ngpy],orderinc);
    linc_a = agrid(ia);
    linc_b = bgrid(ib)';
    linc_c = c(orderinc);
    linc_h = h(orderinc);
    linc_nw = lnw(iab);
end

%% means
Ehours = sum(sum(sum(h.*lmargdist)));
Elabor = sum(sum(sum(llabor.*lmargdist)));
Ewage = sum(sum(sum(lwage.*lmargdist)));
Enetlabinc = sum(sum(sum(lnetlabinc.*lmargdist)));
Egrosslabinc = sum(sum(sum(lgrosslabinc.*lmargdist)));
Enetprofinc = sum(sum(sum(lnetprofinc.*lmargdist)));
Egrossprofinc = sum(sum(sum(lgrossprofinc.*lmargdist)));
Einc = sum(lincgrid.*lincmargdist);
Ea = sum(sum(sum(agrid.*lmargdist)));
Eb = sum(sum(sum(bgrid.*lmargdist)));
Ec = sum(sum(sum(c.*lmargdist)));
Ed = sum(sum(sum(d.*lmargdist)));
Eadjcost = sum(sum(sum(ladjcost.*lmargdist)));
EbP = sum(sum(sum(bgrid(bgrid>0).*lmargdist(:,bgrid>0,:))));
EbN = sum(sum(sum(bgrid(bgrid<0).*lmargdist(:,bgrid<0,:))));
Enw = sum(lnwgrid.*lnwmargdist);

% frac at exactly zero
FRACa0 = lamargdist(1);
if Borrowing; FRACb0 = lbmargdist(ngpbNEG+1); else; FRACb0 = lbmargdist(1); end
FRACnw0 = sum(lnwmargdist(abs(lnwgrid)<1e-8));
FRACb0a0 = sum(lmargdist(abs(la)<1e-8 & abs(lb)<1e-8));
FRACb0aP = sum(lmargdist(la>1e-8 & abs(lb)<1e-8));

% frac zero or close but greater than zero
FRACb0close = interp1(bgrid,lbmargcum,defnbclose*Egrosslabinc);
if Borrowing; FRACb0close = FRACb0close - lbmargcum(ngpbNEG); end
FRACa0close = interp1(agrid,lamargcum,defnaclose*Egrosslabinc);
if ~iteratingtransition && ~calibrating
    FRACnw0close = interp1(lnwgrid,lnwmargcum,(defnbclose+defnaclose)*Egrosslabinc);
    if Borrowing; FRACnw0close = FRACnw0close - interp1(lnwgrid,lnwmargcum,-1e-8); end
end
if ~iteratingtransition % since we don't store gabcum
    FRACb0a0close = interp2(agrid,bgrid,gabcum',defnaclose*Egrosslabinc,defnbclose*Egrosslabinc);
    if Borrowing; FRACb0a0close = FRACb0a0close - interp1(agrid,gabcum(:,ngpbNEG),defnaclose*Egrosslabinc); end
end

FRACbN = sum(lmargdist(lb<-1e-12));
c_b = c.*lmargdist;
Ec_bN = sum(sum(sum(c_b(:,bgrid<-1e-8,:))))/FRACbN;
Ec_b0close = sum(sum(sum(c_b(:,bgrid>-1e-8 & bgrid<defnbclose*Egrosslabinc,:))))/FRACb0close;
Ec_b0far = sum(sum(sum(c_b(:,bgrid>=defnbclose*Egrosslabinc,:))))/(1-FRACb0close-FRACbN);

% percentiles: use cumulative marginal distributions
lpvec = [1 2 5 10 25 50 75 90 95 98 99]/100;
if ~iteratingtransition && (~calibrating || MatchMedianLiq) && (ngpy>1 || deathrate>0)
    % liquid wealth
    iab = diff(lbmargcum) > 0;
    PERCb = interp1(lbmargcum(iab),bgrid(iab),lpvec);
    PERCb(lpvec<=lbmargcum(1)) = bgrid(1);
end
if ~iteratingtransition && (~calibrating || MatchMedianIll || MatchP75Ill) && (ngpy>1 || deathrate>0)
    % illiquid wealth
    iab = diff(lamargcum) > 0;
    PERCa = interp1(lamargcum(iab),agrid(iab),lpvec);
    PERCa(lpvec<=lamargcum(1)) = agrid(1);
end
    
if ~iteratingtransition && ~calibrating && (ngpy>1 || deathrate>0)
    % net worth
    iab = diff(lnwmargcum) > 0;
    PERCnw = interp1(lnwmargcum(iab),lnwgrid(iab),lpvec);
    PERCnw(lpvec<=lnwmargcum(1)) = lnwgrid(1);
    % consumption
    iab = diff(lcmargcum) > 0;
    PERCc = interp1(lcmargcum(iab),lcgrid(iab),lpvec);
    PERCc(lpvec<=lcmargcum(1)) = lcgrid(1);
    % gross income
    iab = diff(lincmargcum) > 0;
    PERCinc = interp1(lincmargcum(iab),lincgrid(iab),lpvec);
    PERCinc(lpvec<=lincmargcum(1)) = lincgrid(1);

    % gini coefficient
    if Ea>0; GINIa = sum(lamargcum.*(1-lamargcum).*adelta) / Ea; else; GINIa = 0; end
    GINIb = sum(lbmargcum.*(1-lbmargcum).*bdelta) / Eb;
    GINInw = sum(lnwmargcum.*(1-lnwmargcum).*lnwdelta) / Enw;
    GINIc = sum(lcmargcum.*(1-lcmargcum).*lcdelta) / Ec;
    GINIinc = sum(lincmargcum.*(1-lincmargcum).*lincdelta) / Einc;

    % in transition use groupings based on steady sate
    
    % statistics conditional on quartile of net worth distributions    
    lnwmargdist = lnwmargdist./lnwdelta;
    [Ea_nwQ,Eb_nwQ,Ec_nwQ,Einc_nwQ] = ConditionalExpectation(nab,lnwgrid,lnwmargdist,lnwdelta,PERCnw(5:7),lnw_a,lnw_b,lnw_c,lnwgrid);
    
    % statistics conditional on quartile of total gross income distributions
    lincmargdist = lincmargdist./lincdelta;
    [Ea_incQ,Eb_incQ,Ec_incQ,Einc_incQ] = ConditionalExpectation(naby,lincgrid,lincmargdist,lincdelta,PERCinc(5:7),linc_a,linc_b,linc_c,lincgrid);
end