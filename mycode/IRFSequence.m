% function IRFSequence(IRFDir)
% CHARACTER    :: lstring*80
% INTEGER    :: it,ipe
% REAL(8) :: lequity

if SolveStickyPriceTransition
    stickytransition = true;
    if OneAssetNoCapital
        IterateTransOneAssetStickyRb
    else
        IterateTransitionStickyRb
    end
    if SaveCumPolicyFnsIRF; CumulativeConsTransition; end
    if ComputeDiscountedMPC; irf = DiscountedMPCTransition(irf,equmINITSS,solnFINALSS); end
    stickytransition = false;
end

OutputFileIRF = [OutputDir,'IRF_',IRFDir,'/NOFS'];
SaveIRFOutput(irf,OutputFileIRF)

%%
if DoPriceExperiments
    for ipe = 1:15
        if SolveStickyPriceTransition
            fprintf('Solving for sticky price transition without ZLB, price experiment %d\n',ipe)
        
            equmTRANS = struct;
            for f = fields(equmINITSS)'
                equmTRANS.(f{1}) = repmat(equmINITSS.(f{1}),1,Ttransition);
            end

            switch ipe
                case 1 % change wage only
                    equmTRANS.wage = irf.equmSTICKY.wage;
                    equmTRANS.netwage = irf.equmSTICKY.netwage;
                    equmTRANS.labtax = irf.equmSTICKY.labtax;

                case 2 % only change profits
                    equmTRANS.profit = irf.equmSTICKY.profit;

                case 3 % only change profits and wage
                    equmTRANS.profit = irf.equmSTICKY.profit;
                    equmTRANS.wage = irf.equmSTICKY.wage;
                    equmTRANS.netwage = irf.equmSTICKY.netwage;
                    equmTRANS.labtax = irf.equmSTICKY.labtax;                 
                
                case 4 % only change rb, rborr only
                    equmTRANS.rb = irf.equmSTICKY.rb;
                    equmTRANS.rborr = irf.equmSTICKY.rborr;

                case 5 % only change ra only
                    equmTRANS.ra = irf.equmSTICKY.ra;

                case 6 % only change illiquid asset drop
                    equmTRANS.illassetdrop = irf.equmSTICKY.illassetdrop;

                case 7 % only change transfers
                    equmTRANS.lumptransfer = irf.equmSTICKY.lumptransfer;

                case 8 % change all 6
                    equmTRANS.wage = irf.equmSTICKY.wage;
                    equmTRANS.netwage = irf.equmSTICKY.netwage;
                    equmTRANS.labtax = irf.equmSTICKY.labtax;
                    equmTRANS.lumptransfer = irf.equmSTICKY.lumptransfer;
                    equmTRANS.rb = irf.equmSTICKY.rb;
                    equmTRANS.rborr = irf.equmSTICKY.rborr;
                    equmTRANS.ra = irf.equmSTICKY.ra;
                    equmTRANS.illassetdrop  = irf.equmSTICKY.illassetdrop;
                    equmTRANS.profit = irf.equmSTICKY.profit;

                 case 9 % change lump transfer by direct effect from government interest payments
                    equmTRANS.lumptransfer = equmINITSS.lumptransfer + (irf.equmSTICKY.rb - equmTRANS.rb)*equmINITSS.govbond;

                 case 10 % change lump transfer, not including direct effect from govt interest payments
                    equmTRANS.lumptransfer = irf.equmSTICKY.lumptransfer - (irf.equmSTICKY.rb - equmTRANS.rb)*equmINITSS.govbond;

                case 11 % change rb, rborr, and change ra by same amount as rb
                    equmTRANS.rb = irf.equmSTICKY.rb;
                    equmTRANS.rborr = irf.equmSTICKY.rborr;
                    equmTRANS.ra = equmINITSS.ra + (irf.equmSTICKY.rb - equmINITSS.rb);

                case 12 % change ra, and change rb, rborr by same amount as ra
                    equmTRANS.ra = irf.equmSTICKY.ra;
                    equmTRANS.rb = irf.equmSTICKY.rb + (irf.equmSTICKY.ra - equmINITSS.ra);
                    equmTRANS.rborr = irf.equmSTICKY.rborr + (irf.equmSTICKY.ra - equmINITSS.ra);

                case 13 % change only proportional labor tax
                    equmTRANS.labtax = irf.equmSTICKY.labtax;
                    equmTRANS.netwage = equmINITSS.wage*(1-irf.equmSTICKY.labtax);

                case 14 % change rb, rborr, and change ra by same amount as rb, and discount eqm profits at implied ra
                    equmTRANS.rb = irf.equmSTICKY.rb;
                    equmTRANS.rborr = irf.equmSTICKY.rborr;
                    equmTRANS.ra = equmINITSS.ra + (irf.equmSTICKY.rb - equmINITSS.rb);
                    if ~DividendFundLumpSum
                        equmTRANS.illassetdrop(:) = 1;
                    elseif DividendFundLumpSum
                        it = Ttransition;
                        lequity = (equmFINALSS.equity + irf.equmSTICKY.dividend(it)*deltatransvec(it)) / (1+deltatransvec(it)*equmTRANS.ra(it));
                        for it = Ttransition-1:-1:1
                            lequity = (lequity + irf.equmSTICKY.dividend(it)*deltatransvec(it)) / (1+deltatransvec(it)*equmTRANS.ra(it));
                        end
                        equmTRANS.illassetdrop(:) = ((1-irf.equmSTICKY.fundlev(1)) * irf.equmSTICKY.capital(1) + lequity) / ((1-equmINITSS.fundlev)*equmINITSS.capital + equmINITSS.equity);
                    end
                    
                case 15 % change rb, rborr, and change ra by same amount as rb, and discount initimake al profits at implied ra
                    equmTRANS.rb = irf.equmSTICKY.rb;
                    equmTRANS.rborr = irf.equmSTICKY.rborr;
                    equmTRANS.ra = equmINITSS.ra + (irf.equmSTICKY.rb - equmINITSS.rb);
                    if ~DividendFundLumpSum
                        equmTRANS.illassetdrop(:) = 1;
                    elseif DividendFundLumpSum
                        it = Ttransition;
                        lequity = (equmFINALSS.equity + equmINITSS.dividend*deltatransvec(it)) / (1+deltatransvec(it)*equmTRANS.ra(it));
                        for it = Ttransition-1:-1:1
                            lequity = (lequity + equmINITSS.dividend*deltatransvec(it)) / (1+deltatransvec(it)*equmTRANS.ra(it));
                        end
                        equmTRANS.illassetdrop(:) = ((1-irf.equmSTICKY.fundlev(1)) * irf.equmSTICKY.capital(1) + lequity) / ((1-equmINITSS.fundlev)*equmINITSS.capital + equmINITSS.equity);
                    end
            end
            Transition

            if SaveCumPolicyFnsIRF; CumulativeConsTransition; end

            irfpriceexp.equmSTICKY = equmTRANS;
            irfpriceexp.statsSTICKY = statsTRANS;
            irfpriceexp.solnSTICKY = solnTRANS;
            irfpriceexp.cumSTICKY = struct('ccum1',ccum1,'ccum2',ccum2,'ccum4',ccum4,...
                                           'dcum1',dcum1,'dcum2',dcum2,'dcum4',dcum4);

            if ComputeDiscountedMPC; irfpriceexp = DiscountedMPCTransition(irfpriceexp,equmINITSS,solnFINALSS); end

        end

        OutputFileIRF = [OutputDir,'IRF_',IRFDir,'/PE',num2str(ipe)];
        SaveIRFOutput(irfpriceexp,OutputFileIRF)
    end    
end



