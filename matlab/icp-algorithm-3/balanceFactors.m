classdef balanceFactors < handle
    %UNTITLED2 This class is used to derive a variety of "balance factors"
    %   Detailed explanation goes here
    
    properties
        img = [];
        env = [];
        snr = [];
        cycles = [];
        snrThresh = .7;
    end
    
    methods
        % constructor to instantiate object with image and envelope
        function obj = balanceFactors(imgIn, envIn)
            if nargin == 0
                error('Must instantiate with image and envelope');
            end
            if length(envIn) ~= size(imgIn, 2)
                error('Envelope length must match image length');
            end
            obj.img = imgIn;
            obj.env = envIn;
        end
        
        function [cycles, SNR] = calculatePeriod(obj)
            % calculates start and end of each good cycle and overall SNR. Returns and stores as local obj
            % returns a vector of the index of the start and end of each systole
            
            % Configurations
            aveConfigs.LP1 = 0.15;
            aveConfigs.aCorrSonoPeaksSpacing = 3;
            aveConfigs.searchWin = 15;
            aveConfigs.corrThresh = .7;
            aveConfigs.S_Dthresh = .25; %.15;
            aveConfigs.fitScoreThresh = 14; %3;
            aveConfigs.slopeThresh = 9.5; %7.5;
            
            % heavily smooth envelope
            [B, A] = butter(4, aveConfigs.LP1, 'low');
            sumSonoSmooth = filtfilt(B, A, obj.env);          
            
            % peak detection
            firstDiff = diff(sumSonoSmooth);
            c = firstDiff(1:end-1).*firstDiff(2:end);
            f = find(c < 0);
            if f(1) == 1
                f(1) = [];
            end
            d =  (firstDiff(f-1) - firstDiff(f)) > 0;
            cycles.sumSonoSmoothPeaks = f(d) + 1;
            cycles.Tvec = diff(cycles.sumSonoSmoothPeaks);
            cycles.T = median(cycles.Tvec);
            aveConfigs.searchWin = round(cycles.T*.15); % x% of period is search window
            
            % reject peaks too close together
            dPeaks = diff(cycles.sumSonoSmoothPeaks);
            f = [find([1000 dPeaks] < round(cycles.T*.4))];
            if ~isempty(f)
                cycles.sumSonoSmoothPeaks(f) = [];    
            end
            % reject peaks outside 15% of median period
            %cycles.goodPeriods = cycles.Tvec <= 1.15*cycles.T & cycles.Tvec >= 0.85*cycles.T;
            %TvecGood = Tvec(goodPeriods);
            %cycles.sumSonoSmoothPeaks = cycles.sumSonoSmoothPeaks(goodPeriods);
                       
            % Look for actual systolic peak (foot) on unsmoothed signal in region of smoothed peak
            %  Room for improvement here
            cycles.systoleStart = [];
            if cycles.sumSonoSmoothPeaks(1) - aveConfigs.searchWin < 1
                cycles.sumSonoSmoothPeaks(1) = [];
            %    cycles.goodPeriods(1) = [];
                cycles.Tvec(1) = [];
            end
            for a = 1:length(cycles.sumSonoSmoothPeaks)
                tmpWin = obj.env(cycles.sumSonoSmoothPeaks(a) - aveConfigs.searchWin:cycles.sumSonoSmoothPeaks(a)); % look back x% of period
                L = length(tmpWin);
                [~, I] = max(tmpWin);
                % get the last value if flat
                if I < L -1
                    while tmpWin(I + 1) == tmpWin(I) && I < L -1
                        I = I + 1;
                    end
                end
                cycles.systolePeak(a) = cycles.sumSonoSmoothPeaks(a) - aveConfigs.searchWin + I - 1;
            end 
            
            cycles.PeriodStart = cycles.systolePeak(1:end-1);
            cycles.PeriodEnd = cycles.systolePeak(2:end);
            
            % check correlations
            
            SNR = 0;
            badEnv = 0;
            removeCycles = [];
            if length(cycles.PeriodStart) > 2
                % if less than 3 cycles then no value in checking and
                % consider this a bad segment
                
                % interpolate each cycle to standard length 
                newLen = 50;
                for a = 1:length(cycles.PeriodStart)
                    tmpEnv = obj.env(cycles.PeriodStart(a) : cycles.PeriodEnd(a));
                    adjEnv(a,:) = interp1(1:length(tmpEnv), tmpEnv, linspace(1, length(tmpEnv), newLen), 'linear','extrap');
                    % median filter
                    mEnv(a,:) = medfilt1(adjEnv(a,:),9);
                    % find intersection points
                    fStart = find(adjEnv(a,:) == mEnv(a,:), 1, 'first');
                    fEnd = find(adjEnv(a,:) == mEnv(a,:), 1, 'last');
                    if isempty(fStart) || isempty(fEnd)
                        removeCycles = [removeCycles; a];
                    else
                        dEnv = abs(adjEnv(a,:) - mEnv(a,:));
                        fitScore = sum(dEnv(fStart:fEnd).^2)/(fEnd - fStart + 1);
                        if fitScore > aveConfigs.fitScoreThresh
                            removeCycles = [removeCycles; a];
                            badEnv = badEnv + 1;
                        end
                    end
                        
                end
                
                SNR = 1 - badEnv/length(cycles.PeriodStart);
                cycles.PeriodStart(removeCycles) = [];
                cycles.PeriodEnd(removeCycles) = [];
                
%                 corrMat = corrcoef(adjEnv');
%                 corrMean = zeros(1,size(corrMat,2));
%                 for b = 1:size(corrMat,2)
%                     tmpCol = corrMat(:,b);
%                     tmpCol(b) = [];
%                     corrMean(b) = mean(tmpCol);
%                 end
%                 fBadCorr = find(corrMean < aveConfigs.corrThresh);
%                 if ~isempty(fBadCorr)
%                     SNR = 1 - length(fBadCorr)/length(corrMean);
%                     cycles.PeriodStart(fBadCorr) = [];
%                     cycles.PeriodEnd(fBadCorr) = [];
%                 else
%                     SNR = 1;
%                 end
            end
            
            % check matches on peak and trough height of remaining cycles
            % also check on length
            vD = []; % diastolic
            vS1 = []; % systolic end
            vL = []; % length (sys - sys)
            if ~isempty(cycles.PeriodStart)
                vS = obj.env(cycles.PeriodStart);
                % get diastolic
                for b = 1:length(vS)
                    peakIdx = cycles.PeriodEnd(b);
                    while obj.env(peakIdx) >= obj.env(peakIdx-1) 
                        peakIdx = peakIdx - 1;
                    end
                    cycles.systoleStart(b) = peakIdx;
                    vD(b) = obj.env(peakIdx);
                    vS1(b) = obj.env(cycles.PeriodEnd(b));
                    vL(b) = cycles.PeriodEnd(b) - cycles.PeriodStart(b) + 1;
                end
                vSmed = median(vS);
                vS1med = median(vS1); % peak at end of peak - peak cycle
                vDmed = median(vD);
                vLmed = median(vL);
                removeCycles = [];
                for b = 1:length(vS)
                    if abs(vS(b) - vSmed)/vSmed > aveConfigs.S_Dthresh
                        removeCycles = [removeCycles; b];
                    elseif abs(vD(b) - vDmed)/vDmed > aveConfigs.S_Dthresh
                        removeCycles = [removeCycles; b];
                    elseif abs(vS1(b) - vS1med)/vS1med > aveConfigs.S_Dthresh
                        removeCycles = [removeCycles; b];
                    elseif abs(vL(b) - vLmed)/vLmed > .2
                        removeCycles = [removeCycles; b];
                    end
                end
                cycles.PeriodStart(removeCycles) = [];
                cycles.PeriodEnd(removeCycles) = [];                
            end
            
            % additional check on sharpness vs roundness of systolic peak
            if ~isempty(cycles.PeriodStart)
                dEnv = diff(obj.env);
                dEnv = [dEnv(1) dEnv];
                removeCycles = [];
                % for each peak at systole start, move backward and ensure
                % sharp slope
                for b = 1:length(cycles.PeriodStart)
                    if abs(dEnv(cycles.PeriodStart(b)-1)) < aveConfigs.slopeThresh && abs(dEnv(cycles.PeriodStart(b)-2)) < aveConfigs.slopeThresh
                        removeCycles = [removeCycles; b];
                    end
                end
                cycles.PeriodStart(removeCycles) = [];
                cycles.PeriodEnd(removeCycles) = [];
            end
            
            % must be at least 3 good cycles to use
            if ~isempty(cycles.PeriodStart)
                if length(cycles.PeriodStart) <= 2
                    SNR = 0;
                    cycles.PeriodStart = [];
                    cycles.PeriodEnd = [];
                end
            else
                SNR = 0;
            end
            obj.snr = SNR;
            obj.cycles = cycles;
        end
        
        function harmonicRatio = getHarmonicRatio(obj, env)
            % calculates the energy in the envelope and obtains ratio of
            % low to high freq energy
            fftLen = 64;
            band1 = 10;
            band2 = 32;
            F = abs(fft(env, fftLen));
            energyBand1 = sum(F(2:band1));
            energyBand2 = sum(F(band1+1:band2));
            
            harmonicRatio = energyBand1/energyBand2;            
        end
        
        function mom = getMoments(obj, img, env, order)
            % method to calculate Nth moment
            %spectralMoment
            [nRows, nCols] = size(img);
            mom = zeros(1, nCols);
            for a = 1:nCols
                yVal = round(env(a));
                if yVal > 1 && yVal <= nRows
                    while isnan(img(yVal, a)) && yVal > 1
                        yVal = yVal - 1;
                    end
                    %mom(a) = trapz( (obj.fAxSigned(1:yVal).^order).*img(1:yVal,a));
                    mom(a) = trapz( ([1:yVal]'.^order).*img(1:yVal,a));
                else
                    mom(a) = NaN;
                end
            end
        end
        
        function balFactors = getBalanceFactors(obj)
            % method to calculate pulsatility index 
            % calculate for each good cycle and then average accross all
            % cycles
            balFactors.snr = 0;
            if obj.snr > obj.snrThresh
                numCycles = length(obj.cycles.PeriodStart);
                for a = 1:numCycles
                    tmpWav = obj.env(obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    tmpImg = obj.img(:,obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    
                    % get PI
                    PI(a) = obj.getPI(tmpWav);
                    
                    % get moments
                    mom1 = obj.getMoments(tmpImg, tmpWav, 1);
                    mom2 = obj.getMoments(tmpImg, tmpWav, 2);
                    
                    % get PI on moments
                    PImom1(a) = obj.getPI(mom1);
                    PImom2(a) = obj.getPI(mom2);
                    
                    % get harmonic ratio on envelope and on moment
                    hrEnv(a) = obj.getHarmonicRatio(tmpWav);
                    hrMom1(a) = obj.getHarmonicRatio(mom1);
                    
                end
                balFactors = obj.getBalFactorsStats('PIenv', balFactors, PI);
                balFactors = obj.getBalFactorsStats('PImom1', balFactors, PImom1);
                balFactors = obj.getBalFactorsStats('PImom2', balFactors, PImom2);
                balFactors = obj.getBalFactorsStats('hrEnv', balFactors, hrEnv);
                balFactors = obj.getBalFactorsStats('hrMom1', balFactors, hrMom1);
                [PI_depth1100, PI_depth2200, PI_depth1111,PI_depth2222,PI_depth1101,PI_depth2201,PI_depth1102,PI_depth2202] = HarmonicsRatio_Envelope(tmpWav,tmpWav);
                [H10ex, H21ex, H32ex, H43ex, H10int, H21int, H32int, H43int, H10ex1, H21ex1, H10int1, H21int1] = HarmonicsRatio_2D_Spectrum(tmpImg, tmpImg);
                
                balFactors.PI_depth1100.mean = PI_depth1100;
                balFactors.PI_depth2200.mean = PI_depth2200;
                balFactors.PI_depth1111.mean = PI_depth1111;
                balFactors.PI_depth2222.mean = PI_depth2222;
                balFactors.PI_depth1101.mean = PI_depth1101;
                balFactors.PI_depth2201.mean = PI_depth2201;
                balFactors.PI_depth1102.mean = PI_depth1102;
                balFactors.PI_depth2202.mean = PI_depth2202;
                
                balFactors.H10ex.mean = H10ex;
                balFactors.H21ex.mean = H21ex;
                balFactors.H32ex.mean = H32ex;
                balFactors.H43ex.mean = H43ex;
                balFactors.H10int.mean = H10int;
                balFactors.H21int.mean = H21int;
                balFactors.H32int.mean = H32int;
                balFactors.H43int.mean = H43int;
                
                balFactors.H10ex1.mean = H10ex1;
                balFactors.H21ex1.mean = H21ex1;
                balFactors.H10int1.mean = H10int1;
                balFactors.H21int1.mean = H21int1;
                
                balFactors.snr = obj.snr;
                
            end
            
        end
        
        function balFactors = getVdStats(obj)
            % method to return the diastolic value, clutter threshold and % diff between Vd and Vclutter 
            % calculate for each good cycle and then average accross all
            % cycles
            balFactors.snr = 0;
            if obj.snr > obj.snrThresh
                numCycles = length(obj.cycles.PeriodStart);
                % get clutter threshold
                meanImg = median(obj.img');
                balFactors.vClutter = find(meanImg > .2*max(meanImg),1,'first')-1;
                if balFactors.vClutter < 1
                    balFactors.vClutter = 1;
                end
                
                for a = 1:numCycles
                    tmpWav = obj.env(obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    tmpImg = obj.img(:,obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    
                    % get Vd
                    Vd(a) = min(tmpWav);
                    
                end
                balFactors = obj.getBalFactorsStats('Vd', balFactors, Vd);
                
                balFactors.snr = obj.snr;
                
            end
            
        end
        
        function balFactors = getPIStats(obj)
            % Similar method to getVdStats but returns the PI for each good cycle and then average accross all
            % cycles
            balFactors.snr = 0;
            if obj.snr > obj.snrThresh
                numCycles = length(obj.cycles.PeriodStart);
                % get clutter threshold
                meanImg = median(obj.img');
                balFactors.vClutter = find(meanImg > .2*max(meanImg),1,'first')-1;
                if balFactors.vClutter < 1
                    balFactors.vClutter = 1;
                end
                
                for a = 1:numCycles
                    tmpWav = obj.env(obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    tmpImg = obj.img(:,obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    
                    % get PI
                    PI(a) = obj.getPI(tmpWav);
                    
                end
                balFactors = obj.getBalFactorsStats('PI', balFactors, PI);
                
                balFactors.snr = obj.snr;
                
            end
            
        end
        
        function balFactors = getNotchStats(obj)
            % Detects and quantifies metrics around a reflected dicrotic
            % notch
            balFactors.snr = 0;
            if obj.snr > obj.snrThresh
                numCycles = length(obj.cycles.PeriodStart);
                obj.cycles.notchTrough = [];
                obj.cycles.notchPeak = [];
                peakRatio = [];
                peakTiming = [];
                               
                for a = 1:numCycles
                    tmpWav = obj.env(obj.cycles.PeriodStart(a):obj.cycles.PeriodEnd(a));
                    % examine diastolic down-slope
                    [minWav, minLoc] = min(tmpWav);
                    obj.cycles.vD(a) = minWav;
                    dSlope = tmpWav(1:minLoc);
                    [pks, loc] = findpeaks(dSlope);
                    if ~isempty(loc)
                        locMins = [];
                        for b = 1:length(loc)
                            locMin = loc(b);
                            while (dSlope(locMin) > dSlope(locMin-1)) && locMin > 1
                                locMin = locMin - 1;
                            end
                            if locMin > 2
                                locMins(b) = locMin;
                            end                                
                        end
                        [~,idx] = max(loc - locMins);
                        obj.cycles.notchPeak(a) = loc(idx) + obj.cycles.PeriodStart(a) - 1;
                        obj.cycles.notchTrough(a) = locMins(idx) + obj.cycles.PeriodStart(a) - 1;
                    
                        peakRatio(a) = (obj.env(obj.cycles.notchPeak(a)) - obj.cycles.vD(a))/(obj.env(obj.cycles.notchTrough(a)) - obj.cycles.vD(a));
                        peakTiming(a) = (obj.cycles.notchPeak(a) - obj.cycles.PeriodStart(a))/(obj.cycles.PeriodEnd(a) - obj.cycles.PeriodStart(a));
                    end
                    peakRatio = peakRatio(peakRatio ~= 0);
                    peakTiming = peakTiming(peakTiming ~= 0);
                        
                end
                balFactors.peakRatio = obj.getBalFactorsStats('peakRatio', balFactors, peakRatio);
                balFactors.peakTiming = obj.getBalFactorsStats('peakTiming', balFactors, peakTiming);
                
                balFactors.snr = obj.snr;
                
            end
            
        end
        
        function balFactors = getBalFactorsStats(obj, factorString, balFactors, data)
            % calculate and append mean, median and std to bal factor data
            % structure
            eval(['balFactors.' factorString '.all = data;' ]);
            eval(['balFactors.' factorString '.mean = mean(data);' ]);
            eval(['balFactors.' factorString '.median = median(data);' ]);
            eval(['balFactors.' factorString '.std = std(data);' ]);
            % obtain mean after removing outliers
            s = std(data);
            f = find(data > 1.96*s + mean(data) | data < mean(data) - 1.96*s);
            if ~isempty(f)
                data(f) = [];
                eval(['balFactors.' factorString '.mean = mean(data);' ]);
            end
                
        end
        
        function PI = getPI(obj, wav)
            % method to calculate pulsatility index 
            L = length(wav);
            PI = (max(wav(1:round(L/3))) - min(wav))/mean(wav);
        end
        
        function getNotch(obj, wav)
            % method to calculate various stats and metrics based on the
            % reflected dicrotic notch and that are less dependent on the
            % peak amplitude which is too variable
            aa = 1;
        end
        
    end
    
end


            