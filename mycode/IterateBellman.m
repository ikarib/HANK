delta = deltass;
	
% set drifts
ltau = 15;
ltau0 = (ra+PerfectAnnuityMarkets*deathrate)*(amax*0.999)^(1-ltau);
adrift = (ra+PerfectAnnuityMarkets*deathrate)*agrid - ltau0*agrid.^ltau;
bdrift = (rborr*(bgrid<=0)+rb*(bgrid>0)+PerfectAnnuityMarkets*deathrate).*bgrid;

if Borrowing && bgrid(1) < -lumptransfer/(rborr+PerfectAnnuityMarkets*deathrate)
	fprintf('Warning: natural borrowing limit violated\n');
end

% Initial Guess
if (~EquilibriumR && ~calibrating) || ...
    (EquilibriumR && neqmiter<=3 && ~calibrating) || ...
    (CalibrateDiscountRate && neqmiter<=3 && ~calibrating) || ...
    (CalibrateRhoAtInitialGuess && neqmiter<=3 && ~calibrating) || ...
    (calibrating && ImposeEqumInCalibration && neqmiter==1 ) || ...
    (calibrating && ~ImposeEqumInCalibration)
	
    switch LaborSupply
        case 0; lh = 1;
        case 1; lh = 1/3;
        case 2; lh = (netwage*ygrid).^frisch;
    end
%     ladrift = (ra + PerfectAnnuityMarkets*deathrate)*agrid;
    lc = (netwage*lh).*ygrid + lumptransfer + (rb+PerfectAnnuityMarkets*deathrate)*bgrid;
    V = utilfn(lc)/(rho+deathrate);
    if LaborSupply==1
        V = V - chi/(rho+deathrate)/(1+1/frisch) * lh.^(1+1/frisch);
    end
    V = repmat(V,ngpa,1,1);
end
ii = 1;
lVdiff = 1.0;
while ii<=maxiter && lVdiff>Vtol
	[Vnew,c,h,s,d,dadj] = HJBUpdate(V,adrift,bdrift,netwage,lumptransfer,labtax,profit,meanlabeff,delta);
	% check for convergence
	lVdiff = max(abs(Vnew(:)-V(:)));
	if Display>1; fprintf('Iteration: %3d    max V change: %.15g\n',ii,lVdiff); end

	% update
	V = Vnew;
	ii = ii+1;
end