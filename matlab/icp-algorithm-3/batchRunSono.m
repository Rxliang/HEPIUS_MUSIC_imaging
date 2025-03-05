% script to batch run all files and generate matching metrics 
pathFolder = sono.getDefaultPath;

load(['pressures_with_ICP_11.9.16.mat']);

d = dir(pathFolder);
isub = [d(:).isdir]; %# returns logical vector
nameFolds = {d(isub).name}';
nameFolds(ismember(nameFolds,{'.','..'})) = [];

numFolders = length(nameFolds);

bAutoTrim = 0; % set to 1 if using auto trim

balFactorsInt = [];
balFactorsExt = [];
balSet = [];

% clear objects
clear balIntObj balExtObj balObj
count = 0;
folderList = [22:64 68:numFolders];
%folderList = 4; %[68:numFolders];

for a = folderList
    numFiles = length(sets{a}.adcFiles);
    oldPressure = -1;
    
    bFirstPressure = 0;
    tmpPressure = [];
    
    pressureIndex = 1;
    withinPressureIndex = 0;
    
    disp(['Folder: ' num2str(a) ' : ' nameFolds{a}]);
    bGoodSet = 0;
    
    icps = [];
    aveIcps = [];
    
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
            
            % auto trim
            if bAutoTrim
                bNegativeEnv = 1;
                trimFactor = 0;
                while bNegativeEnv
                    meanSono = mean(sonoSession.sonoWav.int, 2); % TODO - imrpove auto trim
                    [m, I] = max(meanSono);
                    % make backup
                    sonoImg = trimSono(sonoImg, 0, 'fixed');
                    % trim
                    if I-3 - trimFactor > 0
                        sonoImg = trimSono(sonoImg, I-3 - trimFactor, 'fixed');
                    else
                        bNegativeEnv = 0;
                    end
                    
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
                    
                    if sum(envInt < 1) + sum(envExt < 1) == 0
                        bNegativeEnv = 0;
                    else
                        trimFactor = trimFactor + 1;
                    end
                end
            else
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
            end
            
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
            balIntObj = balanceFactors(sonoImg.int.AR, envInt);
            balIntObj.calculatePeriod();
%             balFactorsInt{pressureIndex, withinPressureIndex} = balIntObj.getBalanceFactors;
            balFactorsInt{pressureIndex, withinPressureIndex} = balIntObj.getVdStats;
            balFactorsInt{pressureIndex, withinPressureIndex}.pressure = oldPressure;
            
            balExtObj = balanceFactors(sonoImg.ext.AR, envExt);
            balExtObj.calculatePeriod();
%             balFactorsExt{pressureIndex, withinPressureIndex} = balExtObj.getBalanceFactors;
            balFactorsExt{pressureIndex, withinPressureIndex} = balExtObj.getVdStats;
            balFactorsExt{pressureIndex, withinPressureIndex}.pressure = oldPressure;

%             % merge and generate icpIndices
%             icps{pressureIndex, withinPressureIndex} = mergeCycles(balIntObj, balExtObj);
%             icps{pressureIndex, withinPressureIndex}.pressure = oldPressure;
%             icps{pressureIndex, withinPressureIndex}.pressure0 = oldPressure0;
%             
%             % factors off averaged sonos
%             [averagedSonoInt, averagedEnvInt] = sonoSession.averageSono(sonoImg.int.AR, envInt, 0);
%             [averagedSonoExt, averagedEnvExt] = sonoSession.averageSono(sonoImg.ext.AR, envExt, 0);
%             
%             aveIcps{pressureIndex, withinPressureIndex} = icpIndices(128 - averagedEnvInt, 128 - averagedEnvExt);
%             aveIcps{pressureIndex, withinPressureIndex}.pressure = oldPressure;
%             aveIcps{pressureIndex, withinPressureIndex}.pressure0 = oldPressure0;
%             balTmp = balanceFactors(sonoImg.ext.AR, envExt);
%             mI(b,:) = balTmp.getMoments(flipud(averagedSonoInt), 128-averagedEnvInt, 1);
%             mE(b,:) = balTmp.getMoments(flipud(averagedSonoExt), 128-averagedEnvExt, 1);
%             pimE(b) = balTmp.getPI(mE(b,:));
%             pimI(b) = balTmp.getPI(mI(b,:));
            
            disp(['Pressure: ' num2str(sonoSession.pressure)])
                       
        end
    end
    
    % update balance stats
    if bGoodSet
        count = count + 1;
        vdStatsInt{count} = scoreVdStats(balFactorsInt);
        vdStatsExt{count} = scoreVdStats(balFactorsExt);
        vdStatsInt{count}.study = nameFolds{a};
        vdStatsExt{count}.study = nameFolds{a};
%         balSetInt{count} = accumulateBalanceFactors(balFactorsInt);
%         balSetExt{count} = accumulateBalanceFactors(balFactorsExt);
%         
%         balObj = matchFactors(balFactorsInt, balFactorsExt);
%         PIenv = balObj.findBalancePressure('PIenv', 0, 1);
%         PImom1 = balObj.findBalancePressure('PImom1', 0, 1);
%         PImom2 = balObj.findBalancePressure('PImom2', 0, 1);
%         hrEnv = balObj.findBalancePressure('hrEnv', 0, 1);
%         hrMom1 = balObj.findBalancePressure('hrMom1', 0, 1);
%         
%         balObj.getBalanceTCD3D;
%         
%         HRE1 = balObj.findBalancePressure('HRE1', 1, 1);
%         HRE2 = balObj.findBalancePressure('HRE2', 1, 0);
%         HRE3 = balObj.findBalancePressure('HRE3', 1, 0);
%         SK0 = balObj.findBalancePressure('SK0', 1, 0);
%         SK1 = balObj.findBalancePressure('SK1', 1, 0);
%         SK2 = balObj.findBalancePressure('SK2', 1, 0);       
%         
%         %balSet.setName = sets{a}.setName;
%         balSet.x(count,:) = [PIenv PImom1 PImom2 hrEnv hrMom1 HRE1 HRE2 HRE3 SK0 SK1 SK2];
%         balSet.y(count) = sets{a}.invICP;
    end
end

save balSet balSet