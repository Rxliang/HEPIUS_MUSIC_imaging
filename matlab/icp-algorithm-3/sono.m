classdef sono < handle
    %sono This class reads in a single .adc file and processes the 10s
    %sonogram for internal and external depths
    %   Detailed explanation goes here
    
    properties
        posPGramFlag=1;
        rawIQ = [];
        sonoWav = [];
        KLD = [];
        env = [];
        binarySono = [];
        binEnv = [];
        binEdge = [];
        vEnv = []; % vittamed envelope already extracted
        fileName = '';
        pressureMap = [];
        pressure = [];
        pressure0 = [];        
        T_s = 10;                
        msWinLen = 150; % 50;
        winOffsetFrac = 0.5; % 0.5;
        HPfilterOrder = 2; %2 HPF order
        cutHPF_Hz = 100; % 100 HPF cutoff
%        startTime_s = 0; 
        % plot properties
        timeAx = [];
        fAxSigned = [];
        fAx=[];
    end
    
    properties (Constant)
        pathConfigFileName = 'pathConfig.xml';
    end
    
    methods 
        
        function obj = loadAdcFile(obj, datFile, dataPath)
            % pop-up dialog to browse to file
            
            % get default folder
            folderName = sono.getDefaultPath;
            % load pressures
            %load([folderName '\'
            %'batchassocpressures_20161101.mat']); % TODO fix to
            %relative path - hardcoded for now 
            sets=[] % jsm mod
            obj.pressureMap = sets;
            
            if nargin < 3
                [datFile dataPath] = uigetfile('*.dat', 'Please browse to .dat file to analyze', folderName);
            end
            filePath = [dataPath datFile];
            obj.rawIQ = load(filePath);
            
            % check to see if correcponding envelope file in directory 
            envFile = [dataPath 'ENV' datFile(4:end)];
            if exist(envFile)
                obj.vEnv.dat = load(envFile);
                obj.vEnv.int = obj.vEnv.dat(:,1);
                obj.vEnv.ext = obj.vEnv.dat(:,2);
                obj.vEnv.t = linspace(0,10,length(obj.vEnv.ext));
            end
            
            obj.fileName = filePath;
            
            % Calculate pressure corresponding to file that has been
            % loaded
            adcFileNum = str2num(datFile(5:7));
            folderNum = [];
            for a = 1:length(sets)
                if ~isempty(strfind(dataPath, sets{a}.setName))
                    folderNum = a;
                end
            end
            if isempty(folderNum)
                %  error('Look into this ...')
                pressure = -2;
                pressure0 = -2;
            else
                
                if ~sets{folderNum}.error
                    f = find(sets{folderNum}.Datafiles <= adcFileNum, 1, 'last');
                    %fHigh = find(sets{folderNum}.Datafiles > adcFileNum, 1, 'first');
                    if ~isempty(f)
                        if mod(f, 2) % odd
                            pressure = sets{folderNum}.pressures((f+1)/2);
                            pressure0 = sets{folderNum}.meanPressure((f+1)/2);
                        elseif sets{folderNum}.Datafiles(f) == adcFileNum
                            pressure = sets{folderNum}.pressures(f/2);
                            pressure0 = sets{folderNum}.meanPressure(f/2);
                        else
                            pressure = NaN;
                            pressure0 = NaN;
                        end
                    else
                        pressure = -100;
                        pressure0 = -100;
                    end
                else
                    pressure = -100;
                    pressure0 = -100;
                end
            end
            obj.pressure = pressure;
            obj.pressure0 = pressure0;
        end

        
        function obj = gensonofwdrevfn(obj, depth, spectrumType, arOrder)
            % depth = 1: process internal depth; depth = 2: process
            % external depth
            % process forward flow only
            
            if nargin < 4
                arOrder = 15;
            end

            I = obj.rawIQ{depth}(:,1);
            Q = obj.rawIQ{depth}(:,2);

            if obj.posPGramFlag
              H1 = imag(hilbert(Q));
              %seq1rPre = I + H1; % reverse
              seq = I - H1; % forward
            else
              seq = I+j*Q;
            end
            
            seqLen = length(seq);
%            fs = seqLen/obj.T_s(depth);
            fs = 1/obj.T_s(depth);

            % remove low frequency clutter
            [B, A] = butter(obj.HPfilterOrder, obj.cutHPF_Hz/(fs/2), 'high');
            %seq1r = filtfilt(B,A,seq1rPre);
            seq = filtfilt(B, A, seq);
            
            preWinLen = obj.msWinLen * fs / 1000;
            winLen = 2^(ceil(log(preWinLen)/log(2)));
            winOffset = max(round(obj.winOffsetFrac * winLen),1); % amount of window

            offset = 1;
            sonoIm = [];

            i = 1;

            while (winLen+offset-1 <= seqLen)
                thisWindow = seq(offset:winLen+offset-1);
                thisWindow = zeropad(thisWindow);
                autocorr = detrend(thisWindow);
                fftLen = length(autocorr);
                switch spectrumType
                    case 'fft'
                        pxx = abs(fft(hamming(fftLen).* autocorr)).^2;
                        posPGram = pxx(1:round(fftLen/2));
                    case 'yule'
                        pxx = pyulear(autocorr, arOrder, 256);
                        posPGram = pxx(1:end-1).^2;
                    case 'burg'
                        pxx = pburg(autocorr, arOrder, 256);
                        posPGram = pxx(1:end-1).^2;
                end
                if obj.posPGramFlag
                  thisLine = flipud(posPGram);
                else
                  thisLine = fftshift(flipud(pxx));
                end
                                    
                sonoIm(:,i) = thisLine;           
                offset = offset + winOffset;
                %                keyboard
                i = i + 1;
            end
            
            winLenOut = length(thisLine);
            % this is the positive freq axis only
            %fAxFull = 0:fs/winLenOut:fs/2*(winLenOut-1)/winLenOut;

            if obj.posPGramFlag
                fAxFull = 0:fs/(winLenOut*2):fs*((winLenOut*2)-1)/(winLenOut*2);
            else
                fAxFull = -fs/2:fs/winLenOut:fs/2-fs/winLenOut;
            end           

            if depth == 1
                obj.sonoWav.int = flipud(sonoIm);
                obj.fAx.int = fAxFull(1:size(obj.sonoWav.int,1));
                obj.timeAx.int = 0:winOffset/fs:(i-2)*winOffset/fs;
                obj.fAxSigned.int = [1:length(obj.fAx.int)]';
            else
                obj.sonoWav.ext = flipud(sonoIm);
                obj.fAx.ext = fAxFull(1:size(obj.sonoWav.ext,1));
                obj.timeAx.ext = 0:winOffset/fs:(i-2)*winOffset/fs;
                obj.fAxSigned.ext = [1:length(obj.fAx.ext)]';                
            end
            
            
            %obj.fAxSigned = [-fliplr(obj.fAx(2:end)) obj.fAx];
            %obj.fAxSigned = [1:length(obj.fAx)]';
            
        end
        
        function obj = statCompMap(obj, depth)
            
            % Automatic Tracing of Blood Flow Velocity in Pulsed Doppler Images
            % Zhe Wang, Student Member IEEE, Greg Slabaugh, Member IEEE Mengchu Zhou, Fellow, IEEE
            % and Tong Fang, Member IEEE
            % calculates the statistical comparison map
            
            if depth == 1
                img = obj.sonoWav.int;
            else
                img = obj.sonoWav.ext;
            end
            
            %img = medfilt2(img, [5 5]);
            img = log10(img);
            % get background
            B = img(end - 10:end, :);
            B = B(:);
            % get image
            I = img(7:17, :);
            I = I(:);
            J = [];
            
            % compute kernel smoothing function estimate for probability function
            [pB, xB] = ksdensity(B);
            [pI, xI] = ksdensity(I);
            pB = pB./(max(pB)*.99);
            pI = pI./(max(pI)*.99);
            
            % calculate divergence map
            for x = 1:size(img, 1)
                for y = 1:size(img, 2)
                    % determine where actual intenstiy is on pdf
                    I1 = img(x, y);
                    fB = find(xB >= I1);
                    if ~isempty(fB)
                        Bi = pB(fB(1));
                    else
                        Bi = pB(end);
                    end
                    fI = find(xI >= I1);
                    if ~isempty(fI)
                        Ii = pI(fI(1));
                    else
                        Ii = pI(end);
                    end
                    
                    J(x,y) = sum(Bi.*log10(Bi./Ii)) + sum(Ii.*log10(Ii./Bi));
                end
            end
            if depth == 1
                obj.KLD.int = J;
            else
                obj.KLD.ext = J;
            end
        end
        
        function env = optimizeGreyThreshold(obj, img, env, greyThresh)
            
            % optimize
            numPeaksOld = peakDetect(env,.3);
            % Move 'slider' to the left in big increments
            shifts = [0.1 0.5];
            numPeaksNew = numPeaksOld;
            N = 2;
            
            while numPeaksNew == numPeaksOld
                greyThresh = 10^(log10(greyThresh) - shifts(2));
                greyThreshOriginal = greyThresh;
                imBin = obj.binarizeSono(img, greyThresh);
                env = obj.getEnvelopeFromBinaryImg(imBin, 1);
                %  env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
                numPeaksNew = peakDetect(env,.3);
                numPeaksOriginal = numPeaksNew;
            end
            
            % Shift one more
            greyThresh = 10^(log10(greyThresh) - shifts(2));
            greyThreshOriginal = greyThresh;
            imBin = obj.binarizeSono(img, greyThresh);
            env = obj.getEnvelopeFromBinaryImg(imBin, 1);
            %  env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
            numPeaksNew = peakDetect(env,.3);
            numPeaksOriginal = numPeaksNew;
            
            % move slider to the right in small increments
            numPeaksOld = numPeaksNew+1;
            numPeaksOld1 = numPeaksOld;
            while numPeaksOld > numPeaksNew || numPeaksOld1 > numPeaksOld || numPeaksNew > 25
                greyThresh = 10^(log10(greyThresh) + shifts(1));
%                 if greyThresh <= greyThreshOriginal
%                     greyThresh = 10^(log10(greyThresh) - N*shifts(2));
%                     N = N + 1;
%                end
                imBin = obj.binarizeSono(img, greyThresh);
                env = obj.getEnvelopeFromBinaryImg(imBin, 1);
                %    env = TVL1denoise(env, figState.fig1.TVRThresh, 100);
                numPeaksOld1 = numPeaksOld;
                numPeaksOld = numPeaksNew;
                numPeaksNew = peakDetect(env,.3);
            end
        end
        
        function bEnv = getEnvelopeFromBinaryImg(obj, img, supressNan)
            
            if nargin < 2
                supressNan = 1;
            end
            
            % simple function to calculate an envelope from a binary image
            nCols = size(img, 2);
            bEnv = NaN(1, nCols);
            for a = 1:nCols
                % find white pixels in each vertical line
                f = find(img(:,a));
                if ~isempty(f)
                    %splitConsecutive = SplitVec(f, 'consecutive', 'split');
                    splitConsecutive = SplitVecSimple(f);
                    lengthConsecutive = cellfun('length', splitConsecutive);
                    severalInRow = find(lengthConsecutive > 2, 1, 'last');
                    if ~isempty(severalInRow)
                        bEnv(a) = cellfun(@(v) v(end), splitConsecutive(severalInRow) );
                    else
                        bEnv(a) = max(f);                            
                    end
                
                end
            end
            
            if supressNan
                % interpolate
                nanEnv = isnan(bEnv);
                t = 1:numel(bEnv);
                bEnv(nanEnv) = interp1(t(~nanEnv), bEnv(~nanEnv), t(nanEnv), 'spline');
                % check if any left over on edges
                f = find(isnan(bEnv));
                if ~isempty(f)
                    if f == 1
                        bEnv(f) = bEnv(2);
                    elseif f == length(bEnv)
                        bEnv(f) = bEnv(f-1);
                    else
                        assert(length(f) > 1, 'Needs fixing');
                    end
                end
            end
            
            % make sure rounded index >=1
            indFix = find(bEnv<0.5);
            bEnv(indFix) = 0.5;
            
        end
        
        function envelope = getEnvelopeGMM(obj, img, fact)
            if nargin < 3
                fact = 1.96;
            end
            
            nCols = size(img, 2);
            envelope = zeros(1, size(img,1));
            
            for a = 1:nCols
                % fit GMM to each col
                gmm =  gmdistribution.fit(img(:,a), 2);
                m1 = gmm.mu(1) + fact*sqrt(gmm.Sigma(:,:,1));
                m2 = gmm.mu(2) + fact*sqrt(gmm.Sigma(:,:,2));
                m = max([m1 m2]);
                f = find(img(:,a) >= m );
                if isempty(f)
                    [maxVal, maxInd] = max(img(:,a));
                    envelope(a) = maxInd;
                else
                    envelope(a) = f(1);
                end
            end
        end
        
        function kalmanEnv = getKalmanEnv(obj, depth, img, z, noiseRatio)
            % uses a Kalman filter to better track the envelope
            % z is the envelope using some other detection 
            
            nCols = size(img, 2);
            img1 = imStretch(img);
            y = z;
            dy = [1 diff(z)];
            %d2y = [1 diff(y)];
            dt = obj.T_s(depth)/nCols;
            
            % initialize
            Xk_prev = [y(1); 0];
            A = [1 dt; 0 1];
            
            nCov = std(dy(2:end))^2;
            Q = nCov;
            
            H = [1 0 ];
            R = Q;
            
            %B = [0.5*dt^2; 1]; % include acceleration term
            B = 0;
            u = 0;
            P = [1  0; 0 1]*nCov;
            
            Xk_buffer = zeros(2,nCols);
            Xk_buffer(:,1) = Xk_prev;
            Z_buffer = zeros(1,nCols);
            
            % run Kalman filter
             for t = 3:nCols-1
                 ps = img1(:,t+1);
                 
                 sig(t) = mean(ps(4:16));
                 noise(t) = mean(ps(end-19:end-4));
                % R = 1/(noise/sig).^noiseRatio;
                 R(t) = (1/(sig/noise))*noiseRatio;
                 
                 Z = y(t+1); 
                 Z_buffer(t+1) = Z;
                 
                 % Kalman iteration
                 P1 = A*P*A' + Q;
                 S = H*P1*H' + R(t);
                 
                 % K is Kalman gain. If K is large, more weight goes to the measurement.
                 % If K is low, more weight goes to the model prediction.
                 K = P1*H'*inv(S);
                 P = P1 - K*H*P1;
                 
                 Xk = A*Xk_prev +B*u + K*(Z-H*A*Xk_prev);
                 Xk_buffer(:,t+1) = Xk;
                 
                 % For the next iteration
                 Xk_prev = Xk;
                 %s(end).R = 1/(noise/sig).^3;%(mean(ps)/thresh)^2;
                 %s(end).z = [z(t); dy(t)]; % create a measurement
                 %s(end+1) = kalmanf(s(end)); % perform a Kalman filter iteration
             end
             kalmanEnv = Xk_buffer;
        end
        
        function mvpEnvelope = getEnvelope1(obj, img)
            
            nCols = size(img, 2);
            snsiEnv = zeros(1, nCols);
            maxFrac = 0.8;
            fMax = size(img, 1);
            x = 0.1;
            
            % get known noise region
            noiseMap = img(end - 15:end, :);
            noiseMap = noiseMap(:);
            thresh = prctile(noiseMap, 85);
            
            for a = 1:nCols
                ps = img(:,a);
                if mean(ps) > thresh
                    ips = cumsum(ps);
                    % strong signal
                    [pmax, mxInd] = max(ps);
                    sigRegionR = mxInd + find(ps(mxInd:end) < pmax*maxFrac, 1) - 1;
                    sigRegionL = find(ps(1:mxInd) < pmax*maxFrac, 1, 'last');
                    linFitSignal = polyfit([sigRegionL sigRegionR], [ips(sigRegionL) ips(sigRegionR)], 1);
                    mSignal = linFitSignal(1);
                    
                    % find knee center
                    % use formula for perpendicular line
                    % http://www.intmath.com/plane-analytic-geometry/perpendicular-distance-point-line.php
                    linFitKneeNoise = polyfit([sigRegionR fMax], [ips(sigRegionR) ips(fMax)], 1);
                    B = 1;
                    A = -1*linFitKneeNoise(1);
                    C = -1*linFitKneeNoise(2);
                    numPts = fMax - sigRegionR + 1;
                    maxDist = 0;
                    maxInd = 1;
                    
                    for b = 1:numPts
                        m = sigRegionR + b -1;
                        n = ips(m);
                        perpDist = abs(A*m + B*n + C)/sqrt(A^2 + B^2);
                        if perpDist > maxDist;
                            maxInd = m;
                            maxDist = perpDist;
                        end
                    end
                    
                    kneeCenter = maxInd;
                    % check to see if its within 1/3 the distance of sigRegionR
                    % to end
                    maxAllowableKnee = ceil((fMax - sigRegionR)/3) + sigRegionR;
                    if kneeCenter > maxAllowableKnee
                        kneeCenter = maxAllowableKnee;
                    end
                    noiseRegion = kneeCenter + (kneeCenter - sigRegionR + 1);
                    
                    linFitNoise = polyfit([noiseRegion fMax], [ips(noiseRegion) ips(fMax)], 1);
                    mNoise = linFitNoise(1);
                    mx = mSignal*x + mNoise*(1 - x);
                    
                    % get slope at all values and find closest match
                    instSlope = [1; diff(ips)];
                    [slopeMatch, indMatch] = min(abs(instSlope(sigRegionR:noiseRegion) - mx));
                    mvpEnvelope(a) = indMatch + sigRegionR - 1;
                else
                    mvpEnvelope(a) = NaN;
                end
            end
            
        end
        
        function imgFilt = simpleFilter(obj, img, filterType, params)
            % function to test some simple filtering ideas
            switch filterType
                case 1
                    % anisotropic filter1
                    imgFilt = exp(-1*(img/params.threshold).^2);
%                     imgFilt = imStretch(imgFilt);
%                     imgFilt = 1 - imgFilt;
                case 2
                    % anisotropic filter2
                    imgFilt = 1./(1 + (img/params.threshold).^2);
                case 3
                    % anisotropic filter3
                    imgFilt = anisodiff2D(img, params.numIters, 1/7, params.threshold, 2);
                case 4
                    % regularization de-noising
                    imgFilt = TVL1denoise(img, params.threshold, 100);
                    
            end
        end
        
        function simpleEnv = getSimpleEnv(obj, img, thresh)
            % uses a simple thresholding on the power spectrum / integrated spectrum
            if nargin < 3
                thresh = 2;
            end
            
            nCols = size(img, 2);
            simpleEnv = zeros(1, nCols);
            
            for a = 1:nCols
                pLine = img(:,a);
                dpsips = smooth(diff(pLine./cumsum(pLine))); % 1st order diff of power spectrum to integrated power spectrum
                noiseThresh = mean(abs(dpsips(end - 20: end-3)));
                f = find( abs(dpsips(1:end - 3)) > noiseThresh*thresh);
                if ~isempty(f)
                    simpleEnv(a) = f(end);
                elseif a > 1
                    simpleEnv(a) = simpleEnv(a - 1);
                else
                    simpleEnv(a) = 0;
                end
            end
        end
        
        function sig = envFilt(obj, img, thresh)
            % apply intelligent total variation filtering to envelope
            % first heavily filter (thresh.low)
            % detect systolic and diastolic and label
            % then moderate filter (thresh.high)
            % apply heavily filtered to diastolic and moderate filter to
            % systolic
            
            if nargin < 3
                thresh.low = .5;
                thresh.high = 1.5;
                thresh.first = 12;
                thresh.second = 25;
            end
            
            env1 = obj.getSimpleEnv(img, thresh.first);   % systolic
            env2 = obj.getSimpleEnv(img, thresh.second);  % diastolic
            
            modFilt = TVL1denoise(env1, thresh.high, 100)*100; % systolic
            heavyFilt = TVL1denoise(env2, thresh.low, 200)*100; % diastolic
            
            % peak detection
            [B, A] = butter(4, .1, 'low');
            sigSmooth = filtfilt(B,A, heavyFilt);
            
            firstDiff = diff(sigSmooth);
            c = firstDiff(1:end-1).*firstDiff(2:end);
            f = find(c < 0);
            if f(1) == 1
                f = f(2:end);
            end
            d = (firstDiff(f-1) - firstDiff(f)) > 0;
            smoothPeaks = f(d) + 1;
            d1 = (firstDiff(f) - firstDiff(f-1)) > 0;
            smoothValleys = f(d1) + 1;
            
            % go back on original image to find actual peak
%             nPeaks = length(smoothPeaks);
%             envPeaks = zeros(1, nPeaks);
%             for a = 1:nPeaks
%                 val = smoothPeaks(a);
%                 sigVal = env1(val);
%                 sigValPre = env1(val-1);
%                 diffVal = sigVal - sigValPre;
%                 while val > 2 && diffVal <= 0
%                     val = val-1;
%                     sigVal = env1(val);
%                     sigValPre = env1(val-1);
%                     diffVal = sigVal - sigValPre;
%                 end
%                 envPeaks(a) = val;
%             end
            
            % replace between smoothed peaks and valleys
            if smoothValleys(1) < smoothPeaks(1)
                smoothValleys = smoothValleys(2:end);
            end
            if smoothPeaks(end) > smoothValleys(end)
                smoothPeaks = smoothPeaks(1:end-1);
            end
            if length(smoothPeaks) ~= length(smoothValleys)
                error('Need to address this');
            end
            sig = modFilt; % systolic
            margin1 = 1;
            margin2 = 5;
            for a = 1:length(smoothValleys)
                if smoothValleys(a)+ margin2 <= length(sig)
                    sig(smoothPeaks(a)-margin1:smoothValleys(a)+margin2) = heavyFilt(smoothPeaks(a)-margin1:smoothValleys(a)+margin2);
                end
            end
            bb = 1;
        end
        
        function snsiEnv = snsiEnvelope(obj, img)
            % This function codes the standard signal noise slope intersection (SNSI) method for
            % envelope detection
            nCols = size(img, 2);
            snsiEnv = zeros(1, nCols);
            maxFrac = 0.7;
            fMax = size(img, 1);
            x = 0.1;
            
            for a = 1:nCols
                ps = img(:,a);
                ips = cumsum(ps);
                % strong signal
                [pmax, mxInd] = max(ps);
                sigRegionR = mxInd + find(ps(mxInd:end) < pmax*maxFrac, 1) - 1;
                sigRegionL = find(ps(1:mxInd) < pmax*maxFrac, 1, 'last');
                linFitSignal = polyfit([sigRegionL sigRegionR], [ips(sigRegionL) ips(sigRegionR)], 1);
                mSignal = linFitSignal(1);
                noiseRegion = round(0.25*(fMax - sigRegionR)) + sigRegionR;  
                linFitNoise = polyfit([noiseRegion fMax], [ips(noiseRegion) ips(fMax)], 1);
                mNoise = linFitNoise(1);
                
                mx = mSignal*x + mNoise*(1 - x);
                %linFitKnee = polyfit([sigRegionR noiseRegion], [ips(sigRegionR) ips(noiseRegion)] ,1);
                %mKnee = linFitKnee(1);
                
                % evaluate x% of the way in 
                %xp = round(x*(noiseRegion - sigRegionR)) + sigRegionR;
                %ipsx = linFitKnee(1)*xp + linFitKnee(2);
                % get slope at all values and find closest match
                instSlope = [1; diff(ips)];
                [slopeMatch, indMatch] = min(abs(instSlope(sigRegionR:noiseRegion) - mx));
                snsiEnv(a) = indMatch + sigRegionR - 1;
                %snsiEnv(a) = find(ips >= ipsx, 1);
                
            end
            
        end
        
        function envSono = createEnvSono(obj, img, env)
            % This function takes a sonogram and a detected envelope and
            % then recreates a sonogram with 0 above the envelope and the
            % original data below the envelope
            
            [nRows, nCols] = size(img);
            envSono = img;
            for a = 1:nCols
                envSono(round(env(a)):nRows,a) = 0;
            end
        end
        
         function cycles = getPeriod(obj, env)
            % returns a vector of the index of the start of each systole
            % for an evenlope series
            
            % Configurations
            aveConfigs.LP1 = 0.15;
            aveConfigs.aCorrSonoPeaksSpacing = 3;
            aveConfigs.searchWin = 15;
            
            [B, A] = butter(4, aveConfigs.LP1, 'low');
            sumSonoSmooth = filtfilt(B, A, env);          
            
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
            aveConfigs.searchWin = round(cycles.T*.35); % x% of period is search window
            
            % reject peaks outside 15% of median period
            cycles.goodPeriods = cycles.Tvec <= 1.15*cycles.T & cycles.Tvec >= 0.85*cycles.T;
            %TvecGood = Tvec(goodPeriods);
            %cycles.sumSonoSmoothPeaks = cycles.sumSonoSmoothPeaks(goodPeriods);
                       
            % Look for start of systolic peak (foot) on unsmoothed signal in region of smoothed peak
            %  Room for improvement here
            cycles.systoleStart = [];
            if cycles.sumSonoSmoothPeaks(1) - aveConfigs.searchWin < 1
                cycles.sumSonoSmoothPeaks(1) = [];
                cycles.goodPeriods(1) = [];
                cycles.Tvec(1) = [];
            end
            for a = 1:length(cycles.sumSonoSmoothPeaks)
                tmpWin = env(cycles.sumSonoSmoothPeaks(a) - aveConfigs.searchWin:cycles.sumSonoSmoothPeaks(a)); % look back x% of period
                [~, I] = min(tmpWin);
                % look forward and make sure sharp systolic peak is found
                if I < length(tmpWin)
                    while tmpWin(I+1) - tmpWin(I) <= 1 & I < length(tmpWin)-1
                        I = I + 1;
                    end
                end
                cycles.systoleStart(a) = cycles.sumSonoSmoothPeaks(a) - aveConfigs.searchWin + I - 1;
            end 
            
        end
        
        function mom = spectralMoment(obj, img, env, order)
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
        
        function [imBin, binThresh] = binarizeSono(obj, img, binThresh)
            % binarizes image. Pass in img from one of the sono class
            % members
            
            % Stretch to convert to gray scale
            im = imStretch(img);
            if nargin <= 2
                %im1 = im(15:end-15,:);
                im1 = im;
                [idx, C] = kmeans(im1(:),2);
                f = find(im1 < min(C));
              %  L1 = length(find(idx == 1));
              %  L2 = length(find(idx == 2));
              %  L = min(L1, L2);
                if std(im1(f)) < 1e-3                    
                    binThresh = prctile(im1(f), 95);
                    if std(im1(f)) < 1e-4
                        binThresh = prctile(im1(f), 85);
                    end
                else
                    binThresh = min(C);
                end
            end
%            im1 = TVL1denoise(im, denoiseThresh, 100);
            
            imBin = im2bw(im, binThresh);
            
        end
        
        function matchStats = matchWaveforms(obj, wave1, wave2, env1, env2, isNormalize)
            % used to generate various matchign stats for waveforms e.g.
            % sum squared error
            if nargin < 5
                isNormalize = 1;
            end
            
            periods1 = obj.getPeriod(env1);
            periods2 = obj.getPeriod(env2);
           
            tol = 3;
            
            % reject cycles
            % find each cardiac cycle to match on
            % align start of systole for env1 and env2 to match
            numGoodPeriods = 0;
            for a = 1:length(periods1.goodPeriods)-1
                if periods1.goodPeriods(a)
                    % identified good period - find match
                    f = find(periods2.systoleStart > periods1.systoleStart(a) - tol & periods2.systoleStart < periods1.systoleStart(a) + tol);
                    if ~isempty(f)
                        if periods2.goodPeriods(f)
                            numGoodPeriods = numGoodPeriods + 1;
                            systoleStart1 = max([periods1.systoleStart(a); periods2.systoleStart(f)]);
                            systoleStart2 = min([periods1.systoleStart(a+1); periods2.systoleStart(f+1)]);
                            
                            
                            cycle1 = wave1(systoleStart1:systoleStart2-1);
                            cycle2 = wave2(systoleStart1:systoleStart2-1);
                            envCycle1 = env1(systoleStart1:systoleStart2-1);
                            envCycle2 = env2(systoleStart1:systoleStart2-1);
                            PI1(numGoodPeriods) = ( max(envCycle1) - min(envCycle1) ) / nanmean(envCycle1);
                            PI2(numGoodPeriods) = ( max(envCycle2) - min(envCycle2) ) / nanmean(envCycle2);
                            %PI1(numGoodPeriods) = ( max(cycle1) - min(cycle1) ) / max(cycle1);
                            %PI2(numGoodPeriods) = ( max(cycle2) - min(cycle2) ) / max(cycle2);
                            if isNormalize
                                cycle1 = cycle1./max(cycle1);
                                cycle2 = cycle2./max(cycle2);
                            end
                            N = length(cycle1);
                            SSE(numGoodPeriods) = nansum( (cycle1 - cycle2).^2 );
                            RMSSD(numGoodPeriods) = 1/N*sqrt( SSE(numGoodPeriods) );
                        end
                    end
                end
            end
            if numGoodPeriods > 0
                matchStats.SSEvec = SSE;
                matchStats.RMSSDvec = RMSSD;
                matchStats.SSE = median(SSE);
                matchStats.RMSSD = median(RMSSD);
                matchStats.PI1 = median(PI1);
                matchStats.PI2 = median(PI2);
                matchStats.numGoodPeriods = numGoodPeriods;
            else
                matchStats.SSEvec = -1;
                matchStats.RMSSDvec = -1;
                matchStats.SSE = -1;
                matchStats.RMSSD = -1;
                matchStats.PI1 = -1;
                matchStats.PI2 = -1;
                matchStats.numGoodPeriods = 0;
            end
        end
        
        function obj = edgeDetect(obj, img, cannyThresh )
            % Canny edge detection on binary image
            if nargin < 3
                cannyThresh = 0.6;
            end
            
            obj.binEdge = edge(img, 'Canny', cannyThresh);
        end
        
        function [averagedSono, averagedEnv] = averageSono(obj, img, env, bDisplay)
            % method to produce averaged sonogram
            % takes a sonogram and envelope as input
            % optional: the period length to use when averaging
            
            if nargin < 4
                bDisplay = 0;
            end
            
            [numRows, numCols] = size(img);
            
            cycles = obj.getPeriod(env);
            T = cycles.T;
            systoleStart = cycles.systoleStart;
            systoleStart(find(cycles.goodPeriods == 0)) = [];
            
            % display
            if bDisplay
                figure;
                imagesc(log10(img)); 
                colormap('jet')
                hold on;    
                plot(numRows - env,'r', 'lineWidth', 2);
                for c = 1:length(systoleStart)
                    line([systoleStart(c) systoleStart(c)],[0 numRows],'Color','g', 'lineWidth', 2)
                end
                pause;
               % close(h);
            end
            
            % set background outside envelope to be 0 before averaging
%             for a = 1:numCols,
%                 img(1: end- env(a)+1, a) = 0;
%             end
            
            % Average pulses.
            % Interpolate to standard length - use length of the average pulse
            newLen = 100;
            adjPulse = zeros(numRows, newLen);
            avPulse = zeros(numRows, newLen);
            avEnv = zeros(1,newLen);
            for a = 1:length(systoleStart)-1
                tmpPulse = img(:, systoleStart(a):systoleStart(a+1));
                tmpEnv = env(systoleStart(a):systoleStart(a+1));
                PulseLen = systoleStart(a+1) - systoleStart(a) + 1;
                
                Y = 1:numRows;
                X = 1:PulseLen;
                Xi = linspace(1, T, newLen);
                for Yi = 1:numRows
                    adjPulse(Yi, :) = interp2(X, Y, tmpPulse, Xi, Yi, 'spline');
                    %adjPulse(Yi, :) = interp1(X, tmpPulse(Yi,:), Xi, 'linear','extrap');
                end
                % interpolate envelope
                adjEnv = interp1(1:PulseLen, tmpEnv, linspace(1, round(T), newLen), 'linear','extrap');
               
                %adjPulse = interp2(1:numRows, 1:PulseLen, tmpPulse, 1:numRows, linspace(1, T, newLen));
                
                %# create interpolant
                %                 [X,Y] = meshgrid(1:pulseLen, 1:numRows);
                %                 F = scatteredInterpolant(X(:), Y(:), tmpPulse(:), 'linear');
                %
                %                 %# interpolate over a finer grid
                %                 [U,V] = meshgrid(linspace(1,T,newLen), 1:numRows);
                %                 adjPulse = F(U,V);
                
                %                 interpolatedPulse = @(x,y) interp2(1:PulseLen,1:numRows,tmpPulse,x,y);
                %                 xIndex = 1;
                %                 adjPulse = zeros(newLen, numRows);
                %                 for xx = linspace(1,T,newLen),
                %                     for yy = 1:numRows,
                %                         adjPulse(xIndex,yy) = interpolatedPulse(xx, yy);
                %                     end
                %                     xIndex = xIndex + 1;
                %                 end
                
                avPulse = avPulse + adjPulse;
                avEnv = avEnv + adjEnv;
            end
            
            averagedSono = abs(avPulse) / a;
            averagedEnv = numRows - (avEnv / a);
            
            % set background outside interpolated envelope to 0
            for b = 1:size(averagedSono, 2),
                averagedSono(1: ceil(averagedEnv(b))+1, b) = NaN;
            end
                       
        end
        
        function plotEnvelopes(obj, env1, env2)
            % function to normalize and plot envelopes
            % env1 = internal
            % env2 = external
            figure;
            h(1) = subplot(3,1,1);
            plot(env1);
            hold on;
            plot(env2, 'r');
            title('Non-normalized envelopes');
            legend('Int', 'Ext');
            
            h(2) = subplot(3,1,2);
            plot(env1./sum(env1));
            hold on;
            plot(env2./sum(env2), 'r');
            title('Normalized envelopes');
            legend('Int', 'Ext');
            
            h(3) = subplot(3,1,3);
            diffEnv = env1./sum(env1) - env2./sum(env2);
            plot(abs(diffEnv));
            title('|Envelope difference|');
            linkaxes(h,'x');
        end
        
        function plotSono(obj, sonoArray)
             % plot multiple sono's time synched on x-axis as an array
             % of structures
             % Pass structures e.g.
             % sonoArray(n).img = obj.sonoWav.int
             % sonoArray(n).isLog = 1
             % sonoArray(n).isGray = 1;
             % sonoArray(n).title = 'Internal Depth'
             % isLog = 1 to plot log plot
             % isGray = 1 to plot gray colormap
             numPlots = size(sonoArray, 2);
             
             for a = 1:numPlots
                 h(a) = subplot(numPlots, 1, a);
                 if sonoArray(a).isLog
                     %imagesc(obj.timeAx, obj.fAxSigned,
                     %log10(sonoArray(a).img))
                     imagesc(obj.timeAx, obj.fAx, log10(sonoArray(a).img))
                     %imagesc(log10(sonoArray(a).img))
                 else
                     imagesc(obj.timeAx, obj.fAx, (sonoArray(a).img))
                     %imagesc((sonoArray(a).img))
                 end
                 
                 % plot the envelope
                 if ~isempty(sonoArray(a).env)
                     hold on
                     envScale = obj.fAx(end) / obj.fAxSigned(end);
                     plot(obj.timeAx, sonoArray(a).env*envScale, 'r', 'lineWidth', 2);
                     %plot(obj.vEnv.t(1:end-9), obj.vEnv.ext(10:end),  'g', 'lineWidth', 2);
                    %plot(sonoArray(a).env, 'r', 'lineWidth', 2);
                     hold off
                 end
                 
                 if sonoArray(a).isGray
                     colormap(gray)
                 end
                 
                 xlabel('t (s)')
                 ylabel('f (Hz)')
                 title(sonoArray(a).title);
                % ylim([-500 1500])
                 axis xy
                 
                 linkaxes(h,'x');
             end
         end
         
              
        function setDefaultPaths(folderName)
            % This function updates the pathConfig.xml file
            % if no file is specified, browse for a file
            if nargin == 0
                folderName = uigetdir('', 'Please browse to the folder with the .adc files');
            end
            
            % create tree
            try
                tree = com.mathworks.xml.XMLUtils.createDocument('pathConfigs');
                clusterNode = tree.getDocumentElement;
                                
                % .dat Data file path
                datFileElement = tree.createElement('dataFilePath');
                datFileElement.appendChild(tree.createTextNode(folderName));
                clusterNode.appendChild(datFileElement);
                                
                xmlFileName = sono.pathConfigFileName;
                xmlwrite(xmlFileName,tree);
            catch
                error(lasterr,' Error');
            end
        end
                
    end
    
    methods (Static)
         function folderName = getDefaultPath
            % This function reads in the pathConfig.xml file to get the
            % default path for browsing
            
            % read in pathConfig.xml to get default path
            % read in and parse xml file
            configFile = parseXML(sono.pathConfigFileName);
            
            % check to see if valid...
            if ~strcmp(configFile.Name, 'pathConfigs')
                error('Invalid xml file');
                return;
            end
            
            for a = 1:length(configFile.Children)
                name = configFile.Children(a).Name;
                switch name
                    case 'dataFilePath'
                        folderName = configFile.Children(a).Children.Data;
                end
            end
            
         end
        
    end
    
end

