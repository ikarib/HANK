% INTEGER     :: it,ii,itfg
% REAL(8)     :: lldK,ldiffB,ldiffK,lminmargcost,lpvgovbc,lpvlumpincr,linitlumpincr
% REAL(8), DIMENSION(Ttransition) :: lcapital,lcapital1,lbond,lfirmdiscount,lrb,lfundbond,lworldbond,lrgov

iteratingtransition = true;

lminmargcost = 0.01;

if Display && stickytransition; fprintf('Solving for sticky price transition\n'); end

% guess capital demand and liquid return

% construct sequence of guesses of capital: assume log linear in capital (constant if a temporary transition)
lldK = log(equmFINALSS.capital/equmINITSS.capital)/Ttransition;
equmTRANS.capital = equmINITSS.capital*exp(lldK*[0 2:Ttransition]);

% construct sequence of guesses of Rb
equmTRANS.pi = repmat(equmINITSS.pi,1,Ttransition);

if IncludeForwardGuideShock
    equmTRANS.rnom = equmINITSS.rnom + [repmat(phifg,1,itfg-1) repmat(phitaylor,1,Ttransition-itfg+1)]*equmTRANS.pi + equmTRANS.mpshock;
else
    equmTRANS.rnom = equmINITSS.rnom + phitaylor*equmTRANS.pi + equmTRANS.mpshock;
end
equmTRANS.rb = equmTRANS.rnom - equmTRANS.pi;

ii = 1;
ldiffK = 1;
ldiffB = 1;
while ii<=maxitertranssticky && max(ldiffK,ldiffB)>toltransition
    % solve for transtion
    tic;Transition;toc
    
    % computed implied equilibrium quantities
    lbond = [statsTRANS.Eb];
    lcapital = ([statsTRANS.Ea] - equmTRANS.equity)./ (1 - equmTRANS.fundlev);
    if ConvergenceRelToOutput
        ldiffK = max(abs(lcapital-equmTRANS.capital)/equmINITSS.output);
        ldiffB = max(abs(lbond-equmTRANS.bond)/equmINITSS.output);
    else
        ldiffK = max(abs(lcapital/equmTRANS.capital - 1));
        ldiffB = max(abs(lbond/equmTRANS.bond - 1));
    end
    if Display
        fprintf('  Transition iter %d:\n',ii)
        fprintf('   K err %.15g,  B err %.15g\n',ldiffK,ldiffB)
        fprintf('   household bond %.15g,  target bond %.15g\n',lbond(2),equmTRANS.bond(2))
    end
    
    % update capital and interest rate
    if ii<maxitertranssticky && max(ldiffK,ldiffB)>toltransition
        equmTRANS.capital(2:Ttransition) = PartialUpdate(Ttransition-1,stepstickytransK,equmTRANS.capital(2:Ttransition),lcapital(2:Ttransition));
        
        lfundbond = -lcapital.*equmTRANS.fundlev;
        lworldbond = -lbond - equmTRANS.govbond - lfundbond;

        lrb = WorldBondInverse2( diff([lworldbond equmFINALSS.worldbond])./(bondadjust*deltatransvec) + lworldbond,equmINITSS.worldbond,equmINITSS.rb,bondelast);
        equmTRANS.rb = PartialUpdate(Ttransition,stepstickytransB,equmTRANS.rb,lrb);
            
    else
        % run distribution stats with full
        iteratingtransition = false;
        Transition
        equmTRANS.capital = lcapital;
        equmTRANS.bond = lbond;
        equmTRANS.rb = lrb;
        
    end
    % inflation
    if IncludeForwardGuideShock
        equmTRANS.pi = (equmTRANS.rb - equmINITSS.rnom - equmTRANS.mpshock) / ([repmat(phifg,1,itfg-1) repmat(phitaylor,1,Ttransition-itfg+1)]-1); % taylor rule
    else
        equmTRANS.pi = (equmTRANS.rb - equmINITSS.rnom - equmTRANS.mpshock) / (phitaylor-1); % taylor rule
    end        
    % nominal interest rates
    equmTRANS.rnom = equmTRANS.rb + equmTRANS.pi; % fisher equn
    % labor
    equmTRANS.labor = [statsTRANS.Elabor];

    save equmTRANS equmTRANS
    save statsTRANS statsTRANS
    ii = ii+1;
end

if stickytransition
    irfpointer.equmSTICKY = equmTRANS;
    irfpointer.statsSTICKY = statsTRANS;
    irfpointer.solnSTICKY = solnTRANS;
end


