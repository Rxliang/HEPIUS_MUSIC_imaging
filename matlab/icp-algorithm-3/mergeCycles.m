function [balFactors] = mergeCycles(balInt, balExt)
%mergeCyles Takes int and ext balObj and finds cycles that are valid in
%both and then generates icpIndices from those valid cycles

% check to see that both have at least one good cycle
if balInt.snr == 0 || balExt.snr == 0
    balFactors.numCycles = 0;
    balFactors.periodStart = [];
    balFactors.periodEnd = [];
    balFactors.icpIndex.index5N = [];
    balFactors.icpIndex.index5 = [];
    balFactors.icpIndex.index2 = [];
    balFactors.icpIndex.index12 = [];
    balFactors.icpIndex.index15 = [];
    balFactors.pressure = [];
    balFactors.pressure0 = [];
    return;   
end

numCyclesInt = length(balInt.cycles.PeriodStart);
numCyclesExt = length(balExt.cycles.PeriodStart);

balFactors.numCycles = min([numCyclesInt numCyclesExt]);

% pick which is smaller
shortStart = balInt.cycles.PeriodStart;
shortEnd = balInt.cycles.PeriodEnd;
longStart = balExt.cycles.PeriodStart;
longEnd = balExt.cycles.PeriodEnd;

if numCyclesExt < numCyclesInt
    shortStart = balExt.cycles.PeriodStart;
    shortEnd = balExt.cycles.PeriodEnd;
    longStart = balInt.cycles.PeriodStart;
    longEnd = balInt.cycles.PeriodEnd;
end

for a = 1:balFactors.numCycles
    % find corresponding
        % biggest start
    start1 = shortStart(a);
    [c, index] = min(abs(longStart - start1));
    start2 = longStart(index);
    periodStart(a) = max([start1 start2]);
    
        % shortest end
    end1 = shortEnd(a);
    [c, index] = min(abs(longEnd - end1));
    end2 = longEnd(index);
    periodEnd(a) = min([end1 end2]);
end

balFactors.periodStart = periodStart;
balFactors.periodEnd = periodEnd;

% calculate icpIndices
for a = 1:balFactors.numCycles
    envInt = balInt.env(periodStart(a):periodEnd(a));
    envExt = balExt.env(periodStart(a):periodEnd(a));
    balFactors.icpIndex(a) = icpIndices(envInt, envExt);
end

balFactors.icpIndex = structofarrays2arrayofstructs(balFactors.icpIndex);

end

