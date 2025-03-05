% script to test variability metrics
bRunAll = 1;

if bRunAll
pathFolder = sono.getDefaultPath;

load(['pressures_with_ICP_11.9.16.mat']);

d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

numFolders = length(nameFolds);

balFactorsInt = [];
balFactorsExt = [];
balSet = [];

% clear objects
clear balIntObj balExtObj balObj
count = 0;

bSelectFolder = 1;
if bSelectFolder
    [selection, ok] = listdlg('PromptString','Select the folder to study:',...
                'SelectionMode','single',...
                'ListString',nameFolds);
    if ok
        folderList = selection;
    else
        return;
    end
else
    folderList = 5;
end

a = folderList;

numFiles = length(sets{a}.adcFiles);
oldPressure = -1;
oldPressure0 = -1;

bFirstPressure = 0;
tmpPressure = [];

pressureIndex = 1;
withinPressureIndex = 0;

disp(['Folder: ' num2str(a) ' : ' nameFolds{a}]);
bGoodSet = 0;
% mI = [];
% mE = [];
% pimI = [];
% pimE = [];
% pres = [];
% eI = [];
% eE = [];
idx = 1;

for b = 1:numFiles
    disp(['File: ' num2str(b) ' : ' sets{a}.adcFiles{b}]);
    % instantiate
    sonoSession = sono;
    sonoSession.msWinLen = 50;
    winOffsetFrac = 0.05;
    
    % load
    sonoSession.loadAdcFile(sets{a}.adcFiles{b}, [pathFolder '\' nameFolds{a} '\Data\']);
    sonoImg = [];
    newPressure = sonoSession.pressure;
    newPressure0 = sonoSession.pressure0;
    
    %valid pressures only
    if ~isnan(sonoSession.pressure) && sonoSession.pressure ~= -100
        bGoodSet = 1;
        %while sonoSession.pressure unchanged
        % calculate spectrum
        % FFT
        %             sonoSession.gensonofwdrevfn(1, 'fft');
        %             sonoSession.gensonofwdrevfn(2, 'fft');
        %             sonoImg.int.fft = sonoSession.sonoWav.int;
        %             sonoImg.ext.fft = sonoSession.sonoWav.ext;
        sonoImg.int.fft = [];
        sonoImg.ext.fft = [];
        % AR
        sonoSession.gensonofwdrevfn(1, 'burg');
        sonoSession.gensonofwdrevfn(2, 'burg');
        sonoImg.int.AR = sonoSession.sonoWav.int;
        sonoImg.ext.AR = sonoSession.sonoWav.ext;
        
        % % Internal envelope detection
        [imBin, binThresh] = sonoSession.binarizeSono(sonoImg.int.AR);
        env1 = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        % optimize
        envInt = sonoSession.optimizeGreyThreshold(sonoImg.int.AR, env1, binThresh);
        
        % External envelope detection
        [imBin, binThresh] = sonoSession.binarizeSono(sonoImg.ext.AR);
        env2 = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
        % optimize
        envExt = sonoSession.optimizeGreyThreshold(sonoImg.ext.AR, env2, binThresh);

        
        if ~bFirstPressure
            bFirstPressure = 1;
            
            oldPressure = newPressure;
            oldPressure0 = newPressure0;
            withinPressureIndex = withinPressureIndex + 1;
        elseif newPressure == oldPressure
            % additional run at same pressure
            
            oldPressure = newPressure;
            oldPressure0 = newPressure0;
            withinPressureIndex = withinPressureIndex + 1;
        else
            % pressure has changed
            oldPressure = newPressure;
            oldPressure0 = newPressure0;
            withinPressureIndex = 1;
            pressureIndex = pressureIndex + 1;
        end
        
        % Generate balance factors
%         balIntObj = balanceFactors(sonoImg.int.AR, envInt);
%         balIntObj.calculatePeriod();
%         balFactorsInt{pressureIndex, withinPressureIndex} = balIntObj.getBalanceFactors;
%         balFactorsInt{pressureIndex, withinPressureIndex}.pressure = oldPressure;
%         
%         balExtObj = balanceFactors(sonoImg.ext.AR, envExt);
%         balExtObj.calculatePeriod();
%         balFactorsExt{pressureIndex, withinPressureIndex} = balExtObj.getBalanceFactors;
%         balFactorsExt{pressureIndex, withinPressureIndex}.pressure = oldPressure;
        
        % factors off averaged sonos
        [averagedSonoInt, averagedEnvInt] = sonoSession.averageSono(sonoImg.int.AR, envInt, 0);
        [averagedSonoExt, averagedEnvExt] = sonoSession.averageSono(sonoImg.ext.AR, envExt, 0);
%         balTmp = balanceFactors(sonoImg.ext.AR, envExt);
%         mI(idx,:) = balTmp.getMoments(flipud(averagedSonoInt), 128-averagedEnvInt, 1);
%         mE(idx,:) = balTmp.getMoments(flipud(averagedSonoExt), 128-averagedEnvExt, 1);
%         pimE(idx) = balTmp.getPI(mE(idx,:));
%         pimI(idx) = balTmp.getPI(mI(idx,:));
%         eI(idx,:) = averagedEnvInt;
%         eE(idx,:) = averagedEnvExt;
        icpIndex(idx) = icpIndices(128-averagedEnvInt, 128-averagedEnvExt);
        
        pres(idx) = newPressure;
        pres0(idx) = newPressure0;
        idx = idx + 1;
        
        disp(['Pressure: ' num2str(sonoSession.pressure)])
        
    end
end

end % if bRunAll

% update balance stats
if bGoodSet
    % re-organize by pressure
    uniquePres = unique(pres);
    uniquePres0 = unique(pres0);
    for x = 1:length(uniquePres)
        f = find(pres == uniquePres(x));
        icpIndex2_tmp = [];
        icpIndex5_tmp = [];
        icpIndex5N_tmp = [];
        icpIndex12_tmp = [];
        icpIndex15_tmp = [];
        if ~isempty(f)
            for b = 1:length(f)
                icpIndex2_tmp(b) = icpIndex(f(b)).index2;
                icpIndex5_tmp(b) = icpIndex(f(b)).index5;
                icpIndex5N_tmp(b) = icpIndex(f(b)).index5N;
                icpIndex12_tmp(b) = icpIndex(f(b)).index12;
                icpIndex15_tmp(b) = icpIndex(f(b)).index15;
            end
        else
            errordlg('Shouldnt get here');
        end
        icpIndex2(x) = median(icpIndex2_tmp);
        icpIndex5(x) = median(icpIndex5_tmp);
        icpIndex5N(x) = median(icpIndex5N_tmp);
        icpIndex12(x) = median(icpIndex12_tmp);
        icpIndex15(x) = median(icpIndex15_tmp);
    end
    
    icps = zeros(15,1);
    Q = zeros(15,1);
    icpFactor = zeros(15,1);
    
    [icps(2), icpFactor(2), ~, Q(2)] = findMinBalance(uniquePres0, uniquePres, icpIndex2, 2);
    [icps(5), icpFactor(5), ~, Q(5)] = findMinBalance(uniquePres0, uniquePres, icpIndex5, 2);
    [icps(6), icpFactor(6), ~, Q(6)] = findMinBalance(uniquePres0, uniquePres, icpIndex5N, 2);
    [icps(12), icpFactor(12), ~, Q(12)] = findMinBalance(uniquePres0, uniquePres, icpIndex12, 2);
    [icps(15), icpFactor(15), ~, Q(15)] = findMinBalance(uniquePres0, uniquePres, icpIndex15, 2);
%     disp('Internal: ');
%     icp = icpAve(mI, pres)
%     disp('External: ');
%     icp = icpAve(mE, pres)
    f = find(Q > 0.8);
    icp.nICP = mean(icps(f));
    icp.invICP = sets{a}.invICP;
    
end
