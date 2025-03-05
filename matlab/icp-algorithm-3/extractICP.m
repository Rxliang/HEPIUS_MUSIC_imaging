% script to run thro xlsx results file and add actual ICP and TCD3D ICP to
% batch structure for ML

% browse for the xlsx file
bLoaded = 1;
if ~bLoaded
    [dataFile dataPath] = uigetfile('*.xls', 'Please browse to .xls trial results file');
    filename = [dataPath dataFile];
    
    % read in file
    [icpData,txt,raw1] = xlsread(filename, 'AARAU ICU Data', 'M3:P98');
    [num,txtData,raw2] = xlsread(filename, 'AARAU ICU Data', 'AK3:AK98');
    
    % load in batch pressure structure
    load('batchassocpressures_20161101.mat');
end

numFiles = size(sets, 2);
for a = 1:numFiles
    % get file name
    fileName = sets{a}.setName;
    % find match
    rowIndex = find(strcmp(fileName, txtData));
    if ~isempty(rowIndex)
        sets{a}.niaICP = icpData(rowIndex,1);
        sets{a}.invICP = icpData(rowIndex,2);
    end
end

save pressures_with_ICP_11.9.16 sets