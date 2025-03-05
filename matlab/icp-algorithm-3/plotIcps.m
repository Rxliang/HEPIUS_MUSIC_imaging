function plotIcps(icpStruct, idx, isAve)

% Takes an icp structure from individual beats and plots on graph
numIcps = size(icpStruct, 1);

if nargin < 3
    isAve = 0;
end

for a = 1:numIcps
    % concatenate
    group1 = cat(2, icpStruct{a,:});
    if ~isAve
        group2 = cat(2, group1.icpIndex);
        fNames = fieldnames(group2);
        data = cat(2, group2.(fNames{idx}));
    else
        fNames = fieldnames(group1);
        data = cat(2, group1.(fNames{idx}));
    end
    plot(icpStruct{a,1}.pressure, data, 'x');
    hold on;
    plot(icpStruct{a,1}.pressure, median(data), 'ro', 'markerSize', 6, 'markerFaceColor', 'k');
    
    pres(a) = icpStruct{a,1}.pressure;
    pres0(a) = icpStruct{a,1}.pressure0;
    
    dataIn(a) = median(data);
end

[ICP, Factor, polynom_data, Q] = findMinBalance(pres, pres0, dataIn, 2);
plot(min(pres):1:max(pres), polynom_data);

title(['ICP = ' num2str(ICP) ' : Q = ' num2str(Q) ' : Factor = ' fNames{idx}])
xlabel('Pressure [mmHg]')