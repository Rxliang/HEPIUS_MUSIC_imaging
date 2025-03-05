% script to organize multi-depth Doppler data into statsInt and statsExt
% structures
rootDir = 'C:\Users\Lance.Myers\Documents\Vittamed\Data\Healthy_volunteer_1_PI_multidepth_analysis\RAW_data';

% get all folders in root dir
folderNamePres = [];
b = 1;
dirs = dir(rootDir);
for a = 3:length(dirs)
    if dirs(a).isdir
        folderNamePres{b} = dirs(a).name;
        b = b+1;
    end
end

numPressures = length(folderNamePres);
% Loop through all pressures
for a = 1:numPressures
    disp(['Pressure: ' folderNamePres{a}]);
    % get depths
    presPath = [rootDir '\' folderNamePres{a}];
    dirs = dir(presPath);
    % assume original dir structure with only one TCD subfolder with data
    for b = 3:length(dirs)
        if dirs(b).isdir
            tcdFolderName{a} = dirs(b).name;
        end
    end
    
    depthPath = [presPath '\' tcdFolderName{a}];
    dirs = dir(depthPath);
    
    c = 1;
    depthFolderName = [];
    for b = 3:length(dirs)
        if dirs(b).isdir
            depthFolderName{c} = dirs(b).name;
            c = c + 1;
        end
    end
    
    % Iterate through each dept
    numDepths = c-1;
    for d = 1:numDepths
        disp(['Depth: ' depthFolderName{d}]);
        % Get number of pressures for given depth
        depthDataPath = [depthPath '\' depthFolderName{d} '\Data'];
        dirData = dir(depthDataPath);
        dataFile = [];
        g = 1;
        for f = 3:length(dirData)
            if strfind(dirData(f).name,'ADC')
                dataFile{g} = dirData(f).name;
                g = g + 1;
            end
        end
        numDataFiles = g-1;
        numDataFiles = numDataFiles - 1; % appears the last file is consistently short and doesnt have sufficient data
        
        % loop through and process each data file
        env1a = [];
        for k = 1:numDataFiles
            disp(['File: ' dataFile{k}]);
            % load file in
            % instantiate
            sonoSession = sono;
            sonoSession.msWinLen = 50;
            winOffsetFrac = 0.05;
            % load
            %currFile = [depthDataPath '\' dataFile{k}];
            sonoSession.loadAdcFile(dataFile{k}, [depthDataPath '\']);
            
            sonoImg = [];
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
            env1a = [env1a env1];
            % optimize
            envInt = sonoSession.optimizeGreyThreshold(sonoImg.int.AR, env1, binThresh);
            
            % External envelope detection
            [imBin, binThresh] = sonoSession.binarizeSono(sonoImg.ext.AR);
            env2 = sonoSession.getEnvelopeFromBinaryImg(imBin, 1);
            % optimize
            envExt = sonoSession.optimizeGreyThreshold(sonoImg.ext.AR, env2, binThresh);
            
            % Use balance factors obj
            balIntObj = balanceFactors(sonoImg.int.AR, envInt);
            balIntObj.calculatePeriod();
            
            % Calculate Vd and PI and store
            %{pressure, depth, datafile} {a,d,k}
            notchStats = balIntObj.getNotchStats;
            if notchStats.snr > 0
                statsInt{a,d,k}.peakRatio = notchStats.peakRatio;
                statsInt{a,d,k}.peakTiming = notchStats.peakTiming;
            end
            statsInt{a,d,k}.Vd = balIntObj.getVdStats;
            statsInt{a,d,k}.PI = balIntObj.getPIStats;
                        
            balExtObj = balanceFactors(sonoImg.ext.AR, envExt);
            balExtObj.calculatePeriod();
            notchStats = balExtObj.getNotchStats;
            if notchStats.snr > 0
                statsExt{a,d,k}.peakRatio = notchStats.peakRatio;
                statsExt{a,d,k}.peakTiming = notchStats.peakTiming;
            end
            statsExt{a,d,k}.Vd = balExtObj.getVdStats;
            statsExt{a,d,k}.PI = balExtObj.getPIStats;
            
            matchObj = matchFactors(balIntObj,balExtObj);
            [wavDist, dPI] = matchObj.compareShapes;
            statsMatch{a,d,k}.wavDist = wavDist;
            statsMatch{a,d,k}.dPI = dPI;
            
        end
        
    end
    
    
end
