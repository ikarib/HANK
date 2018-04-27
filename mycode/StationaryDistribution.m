if (~EquilibriumR && ~calibrating) || ...
    (EquilibriumR && neqmiter<=100 && ~calibrating) || ...
    (CalibrateDiscountRate && neqmiter<=100 && ~calibrating) || ...
    (CalibrateRhoAtInitialGuess && neqmiter<=100 && ~calibrating) || ...
    (calibrating && ImposeEqumInCalibration && neqmiter==1 ) || ...
    (calibrating && ~ImposeEqumInCalibration)

    gmat = zeros(ngpa,ngpb,ngpy);
    ib = 1+~deathrate+Borrowing*ngpbNEG;
    gmat(1,ib,:) = ydist/abdelta(1,ib);
    gmat = reshape(gmat,nab,ngpy);
end

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
B = [1+deltakfe*(sum(AU,2)+deathrate-ymarkovdiag) ... % diagonal value
      -deltakfe*AU.*ABdelta]; % adjust for non-linearly spaced grids

lmat = eye(ngpy) + deltakfe*ymarkovoff';

ldiff = 1.0;
it = 1;
ib = 1+ngpa*ngpbNEG;
while ldiff>KFEtol && it<maxiterKFE
    lgmat = gmat*lmat';
    lgmat(ib,:) = lgmat(ib,:) + deltakfe*deathrate/abdelta(ib)*(abdelta(:)'*gmat);
    % sweep over y
    for iy = 1:ngpy %par
        lgmat(:,iy) = spdiags(B(:,:,iy),[0 -1 1 -ngpa ngpa],nab,nab)\lgmat(:,iy);
    end
   
	ldiff = max(abs(gmat(:)-lgmat(:)));
	gmat = lgmat;
	it = it+1;
	if Display>1; fprintf('KFE iteration: %3d    max g change: %.15g\n',it,ldiff); end
end
gmat(abs(gmat)<1e-50) = 0;
gmat = reshape(gmat,ngpa,ngpb,ngpy);