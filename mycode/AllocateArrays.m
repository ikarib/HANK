fprintf('Allocating arrays for transition\n');

it = 1:Ttransition;

if SolveStickyPriceTransition
    ifs = 1;
    irfstruct(ifs).solnSTICKY(it) = struct;
    if DoPriceExperiments
        irfpriceexp.solnSTICKY(it) = struct;
    end
    if SaveCumPolicyFnsIRF
        irfstruct(ifs).cumSTICKY = struct;
        if DoPriceExperiments
            irfpriceexp.cumSTICKY = struct;
        end
    end
end

fprintf('Finished allocating arrays\n')