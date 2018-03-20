% allocate arrays
AllocateArrays

stickytransition = false;

% set up deltatransvec
deltatransvec = PowerSpacedGrid(Ttransition,deltatransparam,deltatransmin,deltatransmax);
cumdeltatrans = cumsum(deltatransvec);

% equmTRANS = equmINITSS;
equmTRANS.borrwedge = equmINITSS.borrwedge;
equmTRANS.fundlev = repmat(equmINITSS.fundlev,1,Ttransition);
equmTRANS.elast = repmat(equmINITSS.elast,1,Ttransition);
equmTRANS.tfp = repmat(equmINITSS.tfp,1,Ttransition);
equmTRANS.caputil = repmat(equmINITSS.caputil,1,Ttransition);
equmTRANS.labor = repmat(equmINITSS.labor,1,Ttransition);
equmTRANS.labtax = equmINITSS.labtax;
equmTRANS.ra = equmINITSS.ra;

nendtrans = min(50,Ttransition);
irfstruct = struct;

% Monetary policy shock
if IncludeMonetaryShock
    if Display; fprintf('Solving for monetary policy shock IRF\n'); end
    
    equmTRANS.mpshock = equmINITSS.mpshock + MonetaryShockSize * MonetaryShockPers.^[0 cumdeltatrans(1:Ttransition-1)];
    equmTRANS.mpshock(end-nendtrans+1:end) = equmINITSS.mpshock;

    IRFDir = 'Monetary';
    IRFSequence
end

% Forward Guidance shock
if IncludeForwardGuideShock
    if Display; fprintf('Solving for forward guidance shock IRF\n'); end
    
    % forward guidance time
    itfg = find(cumdeltatrans>=ForwardGuideShockQtrs,1);
    equmTRANS.mpshock = equmINITSS.mpshock + ForwardGuideShockSize * MonetaryShockPers.^[zeros(1,itfg) cumsum(deltatransvec(itfg:Ttransition-1))];
    equmTRANS.mpshock(end-nendtrans+1:end) = equmINITSS.mpshock;

    IRFDir = 'ForwardGuide';
    IRFSequence
end