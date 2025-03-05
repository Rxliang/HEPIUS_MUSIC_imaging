% script to comparerevised ICP with new balance factors against invasive
% ICP reference for the Aarau study

load('pressures_with_ICP_11.9.16.mat');
load('newIcps.mat');

for a = 1:newBalance.count
    fileNameBal = newBalance.study{a};
    % make sure it matches to sets{a}
    fileNameSets = sets{a}.setName;
    if isempty(strfind(fileNameBal, fileNameSets))
        % not a match
        errordlg('Check whats going on here')
        return;
    end
    
    % assume matches
    invICP(a) = sets{a}.invICP;
    oldICP(a) = sets{a}.niaICP;
    
    newICP(a) = newBalance.newIcps{a}.ICP;
end