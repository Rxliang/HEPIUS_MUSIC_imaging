function icp = icpAve(icpIn, pres)

% takes in icp and pressure vectors and removes outliers and determines
% variability within and between pressures

icp.P = unique(pres);
numPressures = length(icp.P);

for a = 1:numPressures
    f = find(pres == icp.P(a));
    icp_P = icpIn(f);
    icp.median(a) = median(icp_P);
    icp.std(a) = std(icp_P);
end

icp.var = std(icp.median);
icp.rms = sqrt(sum((icp.median).^2))/length(icp.median);
icp.Mean = mean(icp.median);
icp.stderr1 = std(icp.median)/mean(icp.median);
icp.stderr2 = icp.rms/icp.Mean;
