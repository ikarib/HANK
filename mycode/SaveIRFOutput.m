function SaveIRFOutput(irf,OutputFileIRF)
% Save IRF Output
global SolveStickyPriceTransition SaveTime1PolicyFns SaveCumPolicyFnsIRF
% sticky price transition
if SolveStickyPriceTransition
    if SaveTime1PolicyFns
        for f = fields(irf.solnSTICKY)'
            irf.solnSTICKY.([f{1},'_T1']) = irf.solnSTICKY.(f{1}){1};
            irf.solnSTICKY = rmfield(irf.solnSTICKY,f{1});
        end
        irf.solnSTICKY = rmfield(irf.solnSTICKY,{'AU_T1','gamarg_T1','gbmarg_T1'});
    end    
    if SaveCumPolicyFnsIRF
        for f = fields(irf.cumSTICKY)'
            irf.cumSTICKY.([f{1},'_T1']) = irf.cumSTICKY.(f{1}){1};
            irf.cumSTICKY = rmfield(irf.cumSTICKY,f{1});
        end
    end    
    save(OutputFileIRF,'-struct','irf');
end