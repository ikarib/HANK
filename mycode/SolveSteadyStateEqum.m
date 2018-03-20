% ldiffKN, lKNratio,lKYratio,lwage,lrcapital,lstepKN,ldiffprof,lprofit

fid = fopen(fullfile(OutputDir,'SteadyStateEqumIteration.txt'),'w');

converged = false;
neqmiter = 1;
if calibrating
    lstepKN = 0.2*stepequmss;
else
    lstepKN = stepequmss;
end
ldiffKN = 1;
if ~OneAssetNoCapital
    while (neqmiter<=maxiterequmss && ldiffKN>tolequmss)
load ss neqmiter wage netwage profit ra KNratio KYratio V gmat
        if Display>1
            fprintf('*******************************************\n');
            fprintf(' ITERATION : %.15g\n',neqmiter);
            fprintf(' r guess: %.15g\n',rcapital);
            fprintf('  implied ra: %.15g\n',ra);
            fprintf('  implied borr rate: %.15g\n',rborr);
            fprintf('  implied wage: %.15g\n',wage);
            fprintf('  implied KY ratio: %.15g\n',KYratio);
            fprintf('  implied KN firm: %.15g\n',KNratio);
        end
        
        if initialSS; Grids; end
        IterateBellman
        StationaryDistribution
        DistributionStatistics

        if DividendFundLumpSum
            if DistributeProfitsInProportion
                capital = Ea/(1-fundlev + (1-mc)*(1-corptax)*profdistfrac/(ra*KYratio));
            else
                capital = Ea/(1-fundlev + (1-mc)*(1-corptax)/(ra*KYratio));
            end
        else
            capital = Ea/(1-fundlev);
        end
        
        labor = Elabor;
        equity = Ea - (1-fundlev)*capital;
        lKNratio = capital / labor;
        lKYratio = lKNratio^(1-alpha) / tfp;
        lwage = mc*(1-alpha)*tfp*(lKNratio^alpha);
        lrcapital = mc*alpha/lKYratio;

        ldiffKN = abs(lKNratio/KNratio - 1);
        if Display; fprintf(' Steady state equm iter %2d, K/N error %11.4e\n',neqmiter,lKNratio/KNratio-1); end

        fprintf(fid,'*******************************************');
        fprintf(fid,' ITERATION : %.15g\n',neqmiter);
        fprintf(fid,' r guess: %.15g\n',rcapital);
        fprintf(fid,'  implied ra: %.15g\n',ra);
        fprintf(fid,'  implied borr rate: %.15g\n',rborr);
        fprintf(fid,'  implied wage: %.15g\n',wage);
        fprintf(fid,'  implied KY ratio: %.15g\n',KYratio);
        fprintf(fid,'  actual KY ratio: %.15g\n',lKYratio);
        fprintf(fid,'  implied KN firm: %.15g\n',KNratio);
        fprintf(fid,'  implied KN hh: %.15g\n',lKNratio);
        fprintf(fid,'  relative error: %.15g\n',lKNratio/KNratio-1);

        % update KN ratio
        if neqmiter<=maxiterequmss && ldiffKN>tolequmss
            KNratio = (1-lstepKN)*KNratio+lstepKN*lKNratio;
        else
            KNratio = lKNratio;
        end
        KYratio = (KNratio^(1-alpha)) / tfp;
        profit = (1-mc)*capital/KYratio;
        rcapital = mc*alpha/KYratio;
        wage = mc*(1-alpha)*tfp*(KNratio^alpha);
        netwage = (1-labtax)*wage;
        if DividendFundLumpSum; divrate = 0;
        else; divrate = (1-corptax)*(1-mc)/KYratio; end % outside of steady state include price adjustments
        if DistributeProfitsInProportion; divrate =  profdistfrac*divrate; end
        ra = (rcapital - deprec + divrate - fundlev*rb)/(1-fundlev);
%         if LaborSupply % bug?
%             meanlabeff = (KYratio/KNratio)*3; % scale efficiency units so that average hours would be 1/3 in eqm if cov(z,h)=0
%         else
%             meanlabeff = KYratio^(alpha/(alpha-1)); % scale efficiency units so that output equals 1
%         end

        neqmiter = neqmiter+1;
    end
else

    if DistributeProfitsInProportion
        % iterate on profits
        ldiffprof=1;
        while neqmiter<=maxiterequmss && ldiffprof>tolequmss
            if initialSS; Grids; end
            IterateBellman
            StationaryDistribution
            DistributionStatistics
            labor = Elabor;
            lprofit = (1-mc)*tfp*labor;
            
            ldiffprof = abs(lprofit/profit-1);
            if Display; fprintf(fid,' Steady state equm iter %2d, Profit error %11.4e\n',neqmiter,lprofit/profit-1); end
            
            profit = lprofit;
            neqmiter = neqmiter+1;

        end
    else
        if initialSS; Grids; end
        IterateBellman
        StationaryDistribution
        DistributionStatistics
        labor = Elabor;
        profit = (1-mc)*tfp*labor;
    end
    
    capital = 0;
    equity = 0;
    KYratio = 0;
    KNratio = 0;
end

bond = Eb;
investment = deprec*capital;
priceadjust = 0;
dividend = profit*(1-corptax);
if DistributeProfitsInProportion; dividend = profdistfrac*dividend; end
output = tfp*(capital^alpha)*(labor^(1-alpha));
fundbond = -capital*fundlev;
bondelast = bondelastrelgdp*output;
caputil = 1;
if OneAssetNoCapital; tfpadj = tfp;
else; tfpadj = (tfp^(1+utilelast)*(mc*alpha/rcapital)^(alpha*utilelast))^(1/utilelastalpha); end
taxrev = labtax*wage*labor - lumptransfer + corptax*profit;
if DistributeProfitsInProportion && TaxHHProfitIncome
    taxrev = taxrev + labtax*(1-profdistfrac)*profit*(1-corptax);
end
illassetdrop = 1;

if GovBondResidualZeroWorld
    worldbond = 0;
    govbond = -bond-worldbond-fundbond;
    govexp = taxrev + rb*govbond; 
else
    govbond = -ssdebttogdp*output;
    govexp = taxrev + rb*govbond;
    worldbond = -bond-govbond-fundbond;
end

fclose(fid);