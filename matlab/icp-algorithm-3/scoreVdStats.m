function vdStats = scoreVdStats(vdStruct)

% Function to score various stats around Vd and vClutter for each
% individual cycle

[numPressures, numRepeats] = size(vdStruct);

for a = 1:numPressures
    vdStats.pressure(a) = vdStruct{a,1}.pressure;
    vd = [];
    for b = 1:numRepeats
        if vdStruct{a,b}.snr == 0
            
        else
            vd(b) = vdStruct{a,b}.Vd.mean;
            vClutter(b) = vdStruct{a,b}.vClutter;
        end
        if ~isempty(vd)
            vdStats.vd(a) = mean(vd);
            vdStats.vClutter(a) = mean(vClutter);
            vdStats.diff(a) = (vdStats.vd(a) - vdStats.vClutter(a));
            vdStats.pdiff(a) = vdStats.diff(a)/vdStats.vd(a);            
        else
            vdStats.vd(a) = nan;
            vdStats.vClutter(a) = nan;
            vdStats.diff(a) = nan;
            vdStats.pdiff(a) = nan;
        end
    end
end