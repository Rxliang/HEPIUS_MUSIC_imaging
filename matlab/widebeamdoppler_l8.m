% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use
%
% File name: SetUpGEL3_12DWideBeamDoppler.m - Example of WideBeam doppler
% imaging used in Clnical Application
%
% Description: 
%   Sequence programming file for GEL3-12D Linear array, using WideBeam
%   transmits for 2D (B-mode) and doppler imaging.
%   - this version uses TPC profiles to provide a seperate high voltage
%     level for Doppler. 
%   - Of the 128 transmit channels, only numTx are used, with numTx transmitters on each side of the center element 
%     (where possible). All 128 receive channels are used, although the 
%     element sensitivity cutoff will limit the useful aperture. The receive acquisitions
%     use 200% bandwidth. 
%   - Four Widebeams are used for Doppler acquisition and the aperture for
%     HvMux setting is changed for per Dopper beam
%   - processing is asynchronous with respect to acquisition. 
%   - The image depth, angle, and acquisition PRF can be adjused in Advanced sccript
%   - An outline on the display shows the area undergoing Doppler processing
%   
% Last update:
% 11/28/2016 - test with SW 3.2

clear all

mfile = mfilename;

vtHome = getenv('MUSIC_HOME');
vtData = fullfile(vtHome, 'data');
dataPath = fullfile(vtData, mfile);

if ~exist(dataPath, 'dir')
  mkdir(dataPath)
end

outFilePre = fullfile(dataPath, mfile);
outFileCDIPre = fullfile(dataPath, [mfile '_cdi']);

use312=0;
%Trans.frequency =  10.4167;
Trans.frequency =  demodfreqcheckfn(15.625);

Trans.Bandwidth = [4 15];
demodFreq  =   demodfreqcheckfn(7.8125);

Resource.Parameters.speedOfSound = 1540;    % set speed of sound in m/sec before calling computeTrans
lambda_mm = 1e3*Resource.Parameters.speedOfSound/(Trans.frequency* ...
                                                  1e6);

P.startDepthMm = 0.2;  % startDepth in mm
P.endDepthMm   = 8;  % endDepth in mm
P.txFocusMm    = 4*P.endDepthMm;   % transmit focal pt in
                                   % wavelengths (this is what the
                                   % update in callback does)
P.maxDepthMm   = 10;  % maxDepth for RangeChange and RcvBuffer,
                      % reducing hurts slider UI(9) calc

% Specify system parameters.
Resource.Parameters.numTransmit = 256;          % number of transmit channels.
Resource.Parameters.numRcvChannels = 256;        % number of receive channels.
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = 0;

% Specify Trans structure array.
% elementWidth, ElementPos, and lensCorrection attributes are in
% units specified by Trans.units

if 0
Trans.name = 'GEL8-18i';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans_jon_4p0(Trans);  % computeTrans is used for known transducers.
end

if 1
Trans.name = 'GEL8-18i_64';
Trans.units = 'wavelengths'; % Explicit declaration avoids warning message when selected by default
Trans = computeTrans_jon_4p0(Trans);  % computeTrans is used for known transducers.
end

scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

% Set 2D parameters
P.numTx = round(Trans.numelements/4); % no. of elements in TX aperture.
P.numRays = round(Trans.numelements/4); % no. of widebeam rays in frame
numElementsUsed=Trans.numelements;

P.txFocus = P.txFocusMm*scaleToWvl; 
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;

% Set Doppler parameters
PIndCFI=1; % use same P for all
P.dopNumRays = 6; % no. of rays in frame
P(PIndCFI).dopTxOrgChnl = round([1:P.dopNumRays]*Trans.numelements/(P.dopNumRays+1))
P(PIndCFI).dopDispEle = Trans.numelements; % make CFI full width
P(PIndCFI).dopNumTx = floor(P(PIndCFI).dopTxOrgChnl(1)*2);  % no. of elements in TX aperture

%P.dopDispEle = 128; % middle 128 elements width will be used for Doppler display
P.dopStartDepth = P.startDepth;
P.dopEndDepth = P.endDepth;
P.dopPRIs = 16;
P.dopPRF = 2e+03; % Doppler PRF in Hz.
P.dopAngle = 15*pi/180;
P.pwrThres = 0.14; % orig
P.persfreq = 70;
P.perspwr = 80;
P.cpl = 60;  % level at which 2D overwrites Doppler

m = P.dopNumRays;

% Specify PData structure arrays.
% - 2D PData structure
PData(1).PDelta = [Trans.spacing, 0, 0.5];
% jsm modified width of image
imageWidth_mm = numElementsUsed*Trans.spacingMm; % was 25mm
%sizePre = ceil((imageWidth_mm/lambda_mm)/PData(1).PDelta(3)); %
%Trans.spacingMm);
sizePre = ceil(imageWidth_mm/Trans.spacingMm);
%xBound = -sizePre/2*(lambda_mm*PData(1).PDelta(3));
xBound = -sizePre/2*Trans.spacing;
xAxC = linspace(xBound+Trans.spacing/2, -xBound-Trans.spacing/2, sizePre);
originLam = xAxC(1);

PData(1).Size(1,1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % range, origin and pdelta set PData(1).Size.
%PData(1).Size(1,2) =
%ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(1,2) = sizePre;
PData(1).Size(1,3) = 1;      % single image page
%PData(1).Origin =
%[-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z
%of uppr lft crnr.
PData(1).Origin =  [originLam,0,P.startDepth]; % x,y,z of upper lft crnr.

% Define P.numRays rectangular regions centered on TX beam origins.
% - specify P.numRays default Region structures.
rayDelta = Trans.numelements*Trans.spacing/P.numRays;
% where does this width come from (5*rayDelta)
regionWidth2D = 5*rayDelta; % was  5*rayDelta
PData(1).Region = repmat(struct('Shape',struct('Name','Rectangle','Position',[0,0,0],...
                         'width',...
                         regionWidth2D,'height',P.endDepth)),1,P.numRays);

firstRayLocX = -((Trans.numelements-1)/2)*Trans.spacing;

if 1 | use312
  % this calc is not correct for 6sd truncated beams at edges
  % - set position and width of regions to correspond to transmit origins and beam spacing.
  for i = 1:P.numRays
    PData(1).Region(i).Shape.Position(1) = firstRayLocX + (i-1)*rayDelta;
  end
  PData(1).Region = computeRegions(PData);
end

txAxB = linspace(0, Trans.numelements-1, P.dopNumRays+1);
txAxC = axbound2centre(txAxB);
TxDopOrgX = (txAxC-(Trans.numelements-1)/2)*Trans.spacing;

PData(2) = PData(1);
PData(2).PDelta(1) = Trans.spacing/2;
PData(2).Size(1) = ceil((P.endDepth-P.startDepth)/PData(2).PDelta(3));
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
PData(2).Region = repmat(struct('Shape',struct(...
    'Name','Parallelogram',...
    'Position',[0,0,P.dopStartDepth],... %wil be cahnged later
    'width', P.dopDispEle*Trans.spacing/m,...
    'height', P.dopEndDepth-P.dopStartDepth,...
    'angle',P.dopAngle)),1,m);

for n = 1:m    
    PData(2).Region(n).Shape.Position(1) = TxDopOrgX(n);
end
PData(2).Region = computeRegions(PData(2));

% - PData(3) is used to outline the doppler area
PData(3) = PData(2); PData(3).Region = [];
PData(3).Region.Shape = struct(...
    'Name','Parallelogram',...
    'Position',[PData(3).Origin(1)+PData(3).Size(2)*PData(3).PDelta(1)/2,0,P.dopStartDepth],...
    'width', P.dopDispEle*Trans.spacing,...
    'height', P.dopEndDepth-P.dopStartDepth,...
    'angle',P.dopAngle);
PData(3).Region = computeRegions(PData(3));

% Specify Media object.
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
%maxBufSizePerAcq =
%128*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/128);               
maxBufSizePerAcq = 128*ceil(maxBufLength*8*(Trans.frequency/Trans.frequency)/128);      

% Specify Resources.
% - RcvBuffer(1) is for both 2D and Doppler acquisitions.
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = maxBufSizePerAcq*(P.numRays + m*P.dopPRIs); 
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 20;           % 20 frames allocated for RF acqusitions.
% InterBuffer(1) is for 2D reconstructions.
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;          % one intermediate frame needed for 2D.
% InterBuffer(2) is for Doppler reconstructions.
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;          % one intermediate frame needed for Doppler.
Resource.InterBuffer(2).pagesPerFrame = P.dopPRIs;     % P.dopPRIs pages per ensemble
% ImageBuffer(1) is for 2D image.
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 20;
% ImageBuffer(2) is for Doppler image.
Resource.ImageBuffer(2).datatype = 'double';    % image buffer for Doppler
Resource.ImageBuffer(2).numFrames = 20;
% DisplayWindow is for 2D combined with Doppler
Resource.DisplayWindow(1).Title = 'widebeamdoppler_l8';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)]; % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
Resource.DisplayWindow(1).splitPalette = 1;

% ------Specify structures used in Events------
% Specify Transmit waveform structures.  
% - 2D transmit waveform
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.67,2,1];   % A, B, C, D

%load('ge_l3_12_d_chirp_3p5to10', 'TW');

% - Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TW(2).type = 'parametric';
TW(2).Parameters = [demodFreq,0.67,8,1];   

if strcmp(Trans.units, 'wavelengths')
    scaleToWvl = 1; 
end 

numDelays = Trans.numelements;
apodNotSet =  zeros(1, numDelays);
apodSet =  ones(1, numDelays);

% Specify TX structure array.  
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', apodNotSet, ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Trans.numelements),...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax', 15.0), 1, P.numRays+m);
                             
halfEl = numDelays/2;

% - Set event specific TX attributes.
for n = 1:P.numRays  % specify P.numRays transmit events
    TX(n).Origin = [(firstRayLocX + (n-1)*rayDelta), 0.0, 0.0];
    % Compute available transmit mux aperture
    [Dummy,ce] = min(abs(scaleToWvl*Trans.ElementPos(:,1)-TX(n).Origin(1))); % ce is closest ele to cntr of aper.
    lft = round(ce - halfEl);
    if lft < 1, lft = 1; end
    if lft > numDelays+1, lft = numDelays+1; end
    
    % Compute TX.Apod within mux aperture.
    lft = round(ce - P.numTx/2);
    if lft < 1, lft = 1; end
    rt = round(ce + P.numTx/2);
    if rt > numDelays, rt=numDelays; end
     
    TX(n).Apod(lft:rt) = 1;

    if 0 & ~use312
      % put in center of ray even if truncated
      PData(1).Region(n).Shape.Position(1) = xAxC(round(mean([lft rt])));
    end
    
    % Apply apodization function.
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
end

if 0 & ~use312
  PData(1).Region = computeRegions(PData(1));
end

% -- P.dopNumRays TX structs needed for Doppler
%winNum = 192;
winNum = round(192*Trans.numelements/128);
W = hann(winNum)'; % 128 elememts for hann window

if use312
  aper = [17 49 81 113]; % center 128 elements will be used for
                         % doppler
else
  aper = ones(1,m); % there is only aperture 1 for 6SD
end

ApodRange = (-P.dopNumTx/2+1:P.dopNumTx/2);

for n = 1:m
    TX(P.numRays+n).Origin(1) = TxDopOrgX(n);    
    TX(P.numRays+n).Apod(P.dopTxOrgChnl(n)+ApodRange) = W(winNum/2+ApodRange); 
    TX(P.numRays+n).waveform = 2;
    TX(P.numRays+n).Steer = [P.dopAngle,0.0];
    TX(P.numRays+n).Delay = computeTXDelays(TX(P.numRays+n));
    TX(P.numRays+n).TXPD = computeTXPD(TX(P.numRays+n),PData(2));
end

% Specify TPC structures.
TPC(1).name = '2D';
TPC(1).maxHighVoltage = Trans.maxHighVoltage;
TPC(2).name = 'Doppler';
TPC(2).maxHighVoltage = Trans.maxHighVoltage;

RcvProfile(1).LnaGain = 18; % Profile used for imaging
RcvProfile(2).LnaGain = 24; % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity

% Specify TGC Waveform structure.
% - 2D TGC
%TGC(1).CntrlPts = [0, 271, 498, 617, 767, 903, 1000, 1023];
TGC(1).CntrlPts = [0 856.3112 1023 1023 1023 1023 1023 1023]
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
% - Doppler TGC
%TGC(2).CntrlPts = [0, 271, 498, 617, 767, 903, 1000, 1023];

TGC(2).CntrlPts = [0,272,662,662,662,662,662,662]; % original
TGC(2).CntrlPts = [0 293.6730 714.7482 714.7482 714.7482 714.7482 ...
                   714.7482 714.7482];
TGC(2).CntrlPts = [0 317.0729 771.6993 771.6993 771.6993 771.6993 771.6993 771.6993];
TGC(2).rangeMax = P.dopEndDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

% Specify Receive structure arrays -
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.
maxAcqLength2D = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2)-P.startDepth;
maxAcqLengthDop =  sqrt(P.dopEndDepth^2 + (P.dopDispEle*Trans.spacing)^2)-P.dopStartDepth;
samplesPerWaveDop = 2*TW(2).Parameters(1)/Trans.frequency;

% 30% BW
BPF = ...
[+0.00021 +0.00000 -0.00131 +0.00000 +0.00443 +0.00000 -0.01099 ...
 +0.00000 +0.02228 +0.00000 -0.03876 +0.00000 +0.05957 +0.00000 ...
 -0.08209 +0.00000 +0.10257 +0.00000 -0.11697 +0.00000 +0.12213];

wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
Receive = repmat(struct('Apod', apodSet, ...
                        'startDepth', P.startDepth, ...
                        'endDepth', P.startDepth + wl4sPer128*ceil(maxAcqLength2D/wl4sPer128), ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ... % 200% Bandwidth for 2D
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, (P.numRays+m*P.dopPRIs)*Resource.RcvBuffer(1).numFrames);

% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
    % - Set attributes for each frame.
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays  % acquisitions for 2D
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j; 
    end
    an = 1;                    
    for j = (P.numRays+1):(P.numRays+m*P.dopPRIs)
        % Doppler acquisition
        Receive(j+k).startDepth = P.dopStartDepth;
        Receive(j+k).endDepth = P.dopStartDepth+wl2sPer128*ceil(maxAcqLengthDop/wl2sPer128);
        Receive(j+k).sampleMode = 'BS100BW';
        Receive(j+k).demodFrequency = TW(2).Parameters(1);
        Receive(j+k).InputFilter = BPF;
        Receive(j+k).TGC = 2;
        Receive(j+k).framenum = i;
        Receive(j+k).acqNum = j;        % Doppler acqNums continue after 2D
        if isequal(0, mod((j-P.numRays),P.dopPRIs)), an=an+1; end
    end
end                    
    
% Specify Recon structure arrays.
% - We need two Recon structures, one for 2D, one for Doppler. These will be referenced in the same
%   event, so that they will use the same (most recent) acquisition frame.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', zeros(1,1)), 1, 2);
% - Set Recon values for 2D frame.
Recon(1).RINums(1,1:P.numRays) = (1:P.numRays);  % P.numRays ReconInfos needed for P.numRays angles
k = P.numRays + 1;
% - Set Recon values for Doppler ensemble.
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,1];
Recon(2).ImgBufDest = [0,0];
Recon(2).RINums(1,1:m*P.dopPRIs) = (k:(k+m*P.dopPRIs-1));   %
                                                            % P.dopPRIs ReconInfos needed for Doppler ensemble.

% Define ReconInfo structures.
% - For 2D, we need P.numRays ReconInfo structures for P.numRays steering angles.
% - For Doppler, we need P.dopPRIs ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...    % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum',1, ...
                   'scaleFactor', 0.2, ...
                   'regionnum', 1,...
                   'threadSync',1), 1, P.numRays + m*P.dopPRIs);
% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';

%  - ReconInfos for Doppler ensemble.
ReconInfo(P.numRays+1).Pre = 'clearInterBuf';
for n = 1:m
    for j = 1:P.dopPRIs
        ReconInfo(P.numRays+n+(j-1)*m).mode = 'replaceIQ';
        ReconInfo(P.numRays+n+(j-1)*m).txnum = P.numRays + n;
        ReconInfo(P.numRays+n+(j-1)*m).rcvnum = P.numRays+j+(n-1)*P.dopPRIs;
        ReconInfo(P.numRays+n+(j-1)*m).regionnum = n;
        ReconInfo(P.numRays+n+(j-1)*m).pagenum = j;
        ReconInfo(P.numRays+n+(j-1)*m).scaleFactor = 0.2;
    end
end

% Specify Process structure arrays.
%DopState = 'freq';
DopState = 'power';

Process(1).classname = 'Image';
 Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'pdatanum',1,...    % number of PData structure to use
    'pgain',1.0,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','simple',...
    'persistLevel',20,...
    'interpMethod','4pt',...  %method of interp. (1=4pt)
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','lowerHalf',...
    'display',0,...      % don't display image after processing
    'displayWindow',1};

Process(2).classname = 'Doppler'; % process structure for 1st Doppler ensemble
Process(2).method = 'computeCFIPowerEst';
Process(2).Parameters = {'IntBufSrc',[2,1],...          % number of (Inter?)buffer to process.
    'SrcPages',[3,P.dopPRIs-2],...        % start frame number in source buffer
    'ImgBufDest',[2,-1],...
    'pdatanum',2,...           % number of PData structure
    'prf',P.dopPRF,...           % Doppler PRF in Hz
    'wallFilter','regression',...  % 1 -> quadratic regression
    'pwrThreshold',P.pwrThres,...
    'maxPower',50,...
    'postFilter',1};

Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',2,...   % number of buffer to process.
    'framenum',-1,...   % (-1 => lastFrame)
    'srcData','signedColor',... % type of data to display.
    'pdatanum',2,...    % number of PData structure to use
    'pgain',1.0,...            % pgain is image processing gain
    'reject',2,...      % reject level
    'persistMethod','dynamic',...
    'persistLevel',P.persfreq,...
    'interpMethod','4pt',...  %method of interp. (1=4pt)
    'grainRemoval','none',...
    'processMethod','none',...
    'averageMethod','none',...
    'compressMethod','power',...
    'compressFactor',40,...
    'mappingMethod','upperHalf',...
    'threshold',P.cpl,...
    'display',1,...      % display image after processing
    'displayWindow',1};

% EF1 is external function for ROI plot
Process(4).classname = 'External';
Process(4).method = 'wbdl8_roiplot';
Process(4).Parameters = {'srcbuffer','none'};

procImSaveStart=5;
% save B and 
imgSaveMethods = {'imagesavefn', 'imagecdisavefn'};
for p = procImSaveStart:procImSaveStart+1
  Process(p).classname = 'External';
  Process(p).method = imgSaveMethods{p-procImSaveStart+1};
  Process(p).Parameters = {'srcbuffer','image', ... 
                         'srcbufnum', p-procImSaveStart+1, ...
                         'srcframenum', -1, ...
                         'dstbuffer', 'none'};
end

% Specify SeqControl structure arrays.
% - set receive profile for 2D
SeqControl(1).command = 'setRcvProfile';        
SeqControl(1).argument = 1;  % 
% - time between wide beams 
SeqControl(2).command = 'timeToNextAcq';
SeqControl(2).argument = 200;
% -- Time between 2D acquisition and Doppler ensemble. Set to allow time for profile change.
SeqControl(3).command = 'timeToNextAcq';
SeqControl(3).argument = 2000; % time in usec
% -- Change to TPC Profile 2 (Doppler)
SeqControl(4).command = 'setTPCProfile';
SeqControl(4).condition = 'next';
SeqControl(4).argument = 2;
% - set Receive profile for Doppler
SeqControl(5).command = 'setRcvProfile';
SeqControl(5).argument = 2;  %
% - time between Doppler ensumbles (dopPRF)
SeqControl(6).command = 'timeToNextAcq';
SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
%---------------check Doppler PRF--------------------
tof = ceil(2*Receive(P.numRays+1).endDepth/demodFreq);
if SeqControl(6).argument < tof
    SeqControl(6).argument = tof;
    P.dopPRF = round(1/(tof*1e-06));
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
    fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
    fprintf(['"P.dopPRF" is adjusted to ' num2str(P.dopPRF) '\n']);
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.dopPRF; end
   end
end
%---------------check Doppler PRF--------------------
% -- Change to Profile 1 (2D)
SeqControl(7).command = 'setTPCProfile';
SeqControl(7).condition = 'next';
SeqControl(7).argument = 1;
% - time between frames  
SeqControl(8).command = 'timeToNextAcq';
SeqControl(8).argument = 2000;
% - return to Matlab
SeqControl(9).command = 'returnToMatlab';
% - Jump back to start.
SeqControl(10).command = 'jump';
SeqControl(10).argument = 2;
nsc = 11; % nsc is count of SeqControl objects

% Specify Event structure arrays.
n = 1;

Event(n).info = 'ROI plot';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 4;    % outline doppler region
Event(n).seqControl = 0; % no seqCntrl
n = n+1;

for i = 1:Resource.RcvBuffer(1).numFrames
    rcvInd = (P.numRays + m*P.dopPRIs)*(i-1); % rcvInd keeps track of Receive index increment per frame.
    
    Event(n).info = 'Set RcvProfile for 2D';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 1;
    n = n+1;
    
    % Acquire 2D frame
    for j = 1:P.numRays             
        Event(n).info = 'Acquire 2D Wide Beam';
        Event(n).tx = j;   % use next TX structure.
        Event(n).rcv = rcvInd+j;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = 2;
        n = n+1;
    end
    Event(n-1).seqControl = [3,4];   % replace last 2D acquisition Event's seqControl
    
    Event(n).info = 'Set RcvProfile for Doppler';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 5;
    n = n+1;
    
    % Acquire Doppler ensemble.
    for j = 1:m
        for k = 1:P.dopPRIs
            Event(n).info = 'Acquire Doppler ensemble';
            Event(n).tx = P.numRays+j;      % use next TX structure after 2D.
            Event(n).rcv = rcvInd+P.numRays+k+(j-1)*P.dopPRIs;
            Event(n).recon = 0;      % no reconstruction.
            Event(n).process = 0;    % no processing
            Event(n).seqControl = 6; % TTNA for Doppler
            n = n+1;
        end
    end
    Event(n-1).seqControl = [7,8,nsc]; % replace last Doppler acquisition Event's seqControl
    SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
    nsc = nsc+1;
    
    Event(n).info = 'recons and 2D process';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = [1,2];  % reconstruction for both 2D and Doppler
    Event(n).process = 1;    % process 2D
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;
    
    Event(n).info = 'Doppler processing';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 2;    % IQ to estimate processing
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;
    
    Event(n).info = 'Doppler image display';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 3;    % display processing
    Event(n).seqControl = 0; % no seqCntrl

    
    for p = 1:2
      n=n+1;
      Event(n).info = ['save image buf ' num2str(p)];
      Event(n).tx = 0;
      Event(n).rcv = 0;
      Event(n).recon = 0;
      Event(n).process = [procImSaveStart+p-1];
      Event(n).seqControl = 0;
    end
    
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = 9;
    end
    n = n+1;
    
end

Event(n).info = 'Jump back';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = 10;

% User specified UI Control Elements
% - Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutoffCallback');

% - Doppler Mode Button Group
UI(2).Control = {'UserB4','Style','VsButtonGroup','Title','Doppler Mode','NumButtons',2,'Labels',{'Velocity','Power'}};
UI(2).Callback = text2cell('%DopplerModeCallback');

% - Doppler Power Threshold Slider
UI(3).Control = {'UserB3','Style','VsSlider','Label','DopPwrThres','SliderMinMaxVal',[0.0,1.0,P.pwrThres],...
                 'SliderStep',[0.02,0.1],'ValueFormat','%3.2f'};
UI(3).Callback = text2cell('%DopPowerThreshold');

% - Color Priority Threshold Slider
UI(4).Control = {'UserB2','Style','VsSlider','Label','Color Priority','SliderMinMaxVal',[0,255,P.cpl],...
                 'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%ColorPriorityLevel');

% - Color Persistence Slider
UI(5).Control = {'UserB1','Style','VsSlider','Label','Color Persistence','SliderMinMaxVal',[0,100,P.persfreq],...
                 'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
UI(5).Callback = text2cell('%ColorPersistenceLevel');

% - Replot ROI for doppler
UI(6).Control =  {'UserC4','Style','VsButtonGroup','Title','Doppler Display',...
    'NumButtons',2,'Labels',{'on','off'}};
UI(6).Callback = text2cell('%ReplotROI');

% - PRF adjustment
PRFmin = 500; PRFmax = 4500; stepDiff = PRFmax-PRFmin;
UI(7).Control = {'UserC3','Style','VsSlider','Label','PRF (Hz)','SliderMinMaxVal',[PRFmin,PRFmax,P.dopPRF],...
                 'SliderStep',[50/stepDiff,500/stepDiff],'ValueFormat','%3.0f'};
UI(7).Callback = text2cell('%PRFChange');

% - Steer Angle adjustment
UI(8).Control = {'UserC2','Style','VsSlider','Label','Doppler Angle','SliderMinMaxVal',[-20,20,round(P.dopAngle*180/pi)],...
                 'SliderStep',[1/40,5/40],'ValueFormat','%3.0f'};
UI(8).Callback = text2cell('%SteerAngle');

% - Range Change
MinMaxMm = [5,P.maxDepthMm]; % min max in mm
scaleToWvl = Trans.frequency/Resource.Parameters.speedOfSound*1000;
if isfield(Resource.DisplayWindow(1),'AxesUnits') && strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
    AxesUnit = 'mm';
    MinMaxVal = [MinMaxMm, P.dopEndDepth/scaleToWvl];
else
    AxesUnit = 'wls';
    MinMaxVal = [MinMaxMm*scaleToWvl, P.dopEndDepth];
end

stepBase = MinMaxVal(2)-MinMaxVal(1);
UI(9).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],'SliderMinMaxVal',MinMaxVal,...
                 'SliderStep',[5/stepBase,10/stepBase],'ValueFormat','%3.0f'};
UI(9).Callback = text2cell('%RangeChangeCallback');

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;    

% Save all the structures to a .mat file.
filename = 'MatFiles/widebeamdoppler_l8.mat';
save(filename); 
lslrt(filename)
% VSX
return

% **** Callback routines to be encoded by text2cell function. ****
%SensCutoffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};
assignin('base','Control', Control);
return
%SensCutoffCallback

%DopplerModeCallback - Doppler mode change
P = evalin('base','P');
Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
Process = evalin('base','Process');
Resource = evalin('base','Resource');
hDisplay = Resource.DisplayWindow(1).figureHandle;
currentMap = get(hDisplay,'Colormap');

switch UIState
    case 1  % Velocity mode
        newMap = grayscaleCFImap;
        newMap(1:128,:) = currentMap(1:128,:);
        P.perspwr = get(findobj('Tag','UserB1Slider'),'Value');
        Control(1).Parameters = {'Process',2,'method','computeCFIFreqEst'};
        Control(2).Parameters = {'Process',3,'srcData','signedColor','persistMethod','dynamic','persistLevel',P.persfreq};
        Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
        Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
        set(findobj('tag','UserB1Edit'),'String',num2str(P.persfreq,'%3.0f'));
        set(findobj('tag','UserB1Slider'),'Value',P.persfreq);
        assignin('base','DopState','freq');
        % Set modified Process attributes in base Matlab environment.
        Process(2).method = 'computeCFIFreqEst';
        for k = 1:2:length(Process(3).Parameters)
            if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'signedColor';
            elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'dynamic';
            elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = P.persfreq;
            end
        end
    case 2  % Power mode
        newMap = grayscaleCPAmap;
        newMap(1:128,:) = currentMap(1:128,:);
        for k = 1:2:length(Process(3).Parameters)
            if strcmp(Process(3).Parameters{k},'persistLevel'), P.persfreq = Process(3).Parameters{k+1}; end
        end
        Control(1).Parameters = {'Process',2,'method','computeCFIPowerEst'};
        Control(2).Parameters = {'Process',3,'srcData','unsignedColor','persistMethod','simple','persistLevel',P.perspwr};
        Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
        Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
        set(findobj('tag','UserB1Edit'),'String',num2str(P.perspwr,'%3.0f'));
        set(findobj('tag','UserB1Slider'),'Value',P.perspwr);
        assignin('base','DopState','power');
        Process(2).method = 'computeCFIPowerEst';
        for k = 1:2:length(Process(3).Parameters)
            if strcmp(Process(3).Parameters{k},'srcData'), Process(3).Parameters{k+1} = 'unsignedColor';
            elseif strcmp(Process(3).Parameters{k},'persistMethod'), Process(3).Parameters{k+1} = 'simple';
            elseif strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = P.perspwr;
            end
        end
end
assignin('base','P',P);
assignin('base','Process',Process);
assignin('base','Control', Control);

% If PTool window is open, adjust all uicontrols
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    posPTool = get(hPTool,'position');
    PTool;
    set(findobj('tag','ProcessTool'),'position',posPTool);
end
return
%DopplerModeCallback

%DopPowerThreshold - Doppler Power change
Process = evalin('base','Process');
for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'pwrThreshold'), Process(2).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',2,'pwrThreshold',UIValue};
assignin('base','Control', Control);
%DopPowerThreshold

%ColorPriorityLevel - Color Priority change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(3).Parameters)
    if strcmp(Process(3).Parameters{k},'threshold'), Process(3).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',3,'threshold',UIValue};
assignin('base','Control', Control);
%ColorPriorityLevel

%ColorPersistenceLevel - Color Persistence change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
for k = 1:2:length(Process(3).Parameters)
    if strcmp(Process(3).Parameters{k},'persistLevel'), Process(3).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.persistLevel.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',3,'persistLevel',UIValue};
assignin('base','Control', Control);

% If PTool window is open, adjust persistLevel1 in Process(3)
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    hPNum = findobj('tag','processNum');
    if isequal(get(findobj('tag','processNum'),'Value'),3)
        set(findobj('tag','persistSlider1'),'Value',UIValue);
        set(findobj('tag','persistValue1'),'String',num2str(UIValue));
    end
end
return
%ColorPersistenceLevel

%ReplotROI
Process = evalin('base','Process');
if UIState == 1
    evalin('base','set(ROIHandle,''Visible'',''on'');');
    for k = 1:3:length(Process(1).Parameters)
        if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 0; end
    end
    for k = 1:3:length(Process(3).Parameters)
        if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 1; end
    end
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control(1).Command = 'set&Run';
    Control(2).Command = 'set&Run';
    Control(1).Parameters = {'Process',1,'display',0};
    Control(2).Parameters = {'Process',3,'display',1};
else
    evalin('base','set(ROIHandle,''Visible'',''off'');');
    for k = 1:3:length(Process(1).Parameters)
        if strcmp(Process(1).Parameters{k},'display'), Process(1).Parameters{k+1} = 1; end
    end
    for k = 1:3:length(Process(3).Parameters)
        if strcmp(Process(3).Parameters{k},'display'), Process(3).Parameters{k+1} = 0; end
    end
    % Set Control.Command to set Image.persistLevel.
    Control = evalin('base','Control');
    Control(1).Command = 'set&Run';
    Control(2).Command = 'set&Run';
    Control(1).Parameters = {'Process',1,'display',1};
    Control(2).Parameters = {'Process',3,'display',0};
end
assignin('base','Control', Control);
return
%ReplotROI

%PRFChange - PRF
P = evalin('base','P');
P.dopPRF = UIValue;
Trans = evalin('base','Trans');
Process = evalin('base','Process');
demodFreq = evalin('base','demodFreq');
SeqControl = evalin('base','SeqControl');
m = P.dopNumRays;
SeqControl(6).argument = round(1/(P.dopPRF*1e-06));

%---------------check Doppler PRF--------------------
currentDepth = evalin('base','Receive(P.numRays+1).endDepth');
tof = ceil(2*currentDepth/demodFreq);
if SeqControl(6).argument < tof
    SeqControl(6).argument = tof;
    P.dopPRF = round(1/(tof*1e-06));
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
    fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
    fprintf(['"dopPRF" is adjusted to ' num2str(P.dopPRF) '\n']);
    UI = evalin('base','UI');
    set(UI(7).handle(2),'Value',P.dopPRF);
    set(UI(7).handle(3),'String',num2str(P.dopPRF));
end
%---------------check Doppler PRF--------------------

for k = 1:2:length(Process(2).Parameters)
    if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.dopPRF; end
end

assignin('base','P',P);
assignin('base','Process',Process);
assignin('base','SeqControl',SeqControl);

Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',2,'prf',P.dopPRF};
Control(2).Command = 'update&Run';
Control(2).Parameters = {'SeqControl'};
assignin('base','Control',Control);

return
%PRFChange

%SteerAngle - Dopler Angle will change TX beam and PData
P = evalin('base','P');
TX = evalin('base','TX');
PData = evalin('base','PData');

P.dopAngle = UIValue * pi/180;
assignin('base','P',P);

for n = 1:P.dopNumRays
    PData(2).Region(n).Shape.angle = P.dopAngle;
end
PData(3).Region.Shape.angle = P.dopAngle;
PData(2).Region = computeRegions(PData(2));
PData(3).Region = computeRegions(PData(3));
assignin('base','PData',PData);

for n = 1:P.dopNumRays    
    TX(P.numRays+n).Steer = [P.dopAngle,0.0];
    TX(P.numRays+n).Delay = computeTXDelays(TX(P.numRays+n));
    TX(P.numRays+n).TXPD = computeTXPD(TX(P.numRays+n),PData(2));
end
assignin('base','TX',TX);

Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','Recon','TX'};
assignin('base','Control',Control);

return
%SteerAngle

%RangeChangeCallback - range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end

% P.endDepth
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P = evalin('base','P');
P.endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P.endDepth = UIValue*scaleToWvl;
    end
end
P.dopEndDepth = P.endDepth;
P.txFocus = 4*P.dopEndDepth;

% PData
PData = evalin('base','PData');
m = P.dopNumRays;

PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % rows
PData(2) = PData(1); % default is the same as PData(1), but the user may want to change the PDelta
% - Doppler PData structure for TXPD calculation
PData(2) = PData(1);
PData(2).PDelta(1) = Trans.spacing/2;
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1)); % cols
PData(2).Region = repmat(struct('Shape',struct(...
    'Name','Parallelogram',...
    'Position',[0,0,P.dopStartDepth],... %wil be cahnged later
    'width', P.dopDispEle*Trans.spacing/m,...
    'height', P.dopEndDepth-P.dopStartDepth,...
    'angle',P.dopAngle)),1,m);

TxDopOrgX = evalin('base','TxDopOrgX');
for n = 1:m    
    PData(2).Region(n).Shape.Position(1) = TxDopOrgX(n);
end

% - PData(3) is used for doppler processing and outline the doppler area
PData(3) = PData(2); PData(3).Region = [];
PData(3).Region.Shape = struct(...
    'Name','Parallelogram',...
    'Position',[PData(3).Origin(1)+PData(3).Size(2)*PData(3).PDelta(1)/2,0,P.dopStartDepth],...
    'width', P.dopDispEle*Trans.spacing,...
    'height', P.dopEndDepth-P.dopStartDepth,...
    'angle',P.dopAngle);
for n = 1:3
    PData(n).Region = computeRegions(PData(n));
end
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

% Update TXPD data of TX structures.
TX = evalin('base','TX');
h = waitbar(0,'Program TX parameters, please wait!');

for i = 1:P.numRays
    TX(i).focus = P.txFocus;
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(1));
    waitbar(i/(P.numRays+m))
end
for i = P.numRays+1:P.numRays+m
    TX(i).focus = P.txFocus;
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(2));
    waitbar(i/(P.numRays+m))
end
close(h)
assignin('base','TX',TX);

% Receive
Receive = evalin('base', 'Receive');
demodFreq = evalin('base','demodFreq');

maxAcqLngth2D = sqrt(P.endDepth^2 + (Trans.numelements*Trans.spacing)^2) - P.startDepth;
maxAcqLngthDop =  sqrt(P.dopEndDepth^2 + (P.dopDispEle*Trans.spacing)^2) - P.dopStartDepth;
samplesPerWaveDop = 2*demodFreq/Trans.frequency;

wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4 smpls per wave round trip.
wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128 sample block for 2 smpls per wave round trip.
for i = 1:Resource.RcvBuffer(1).numFrames
    k = (P.numRays + m*P.dopPRIs)*(i-1); % k keeps track of Receive index increment per frame.
    for j = 1:P.numRays  % acquisitions for 2D        
        Receive(j+k).endDepth = P.startDepth + wl4sPer128*ceil(maxAcqLngth2D/wl4sPer128); 
    end
    for j = (P.numRays+1):(P.numRays+P.dopPRIs*m); % Doppler acquisition        
        Receive(j+k).endDepth = P.dopStartDepth + wl2sPer128*ceil(maxAcqLngthDop/wl2sPer128);        
    end
end
assignin('base','Receive',Receive);

% PRF
%---------------check Doppler PRF--------------------
Process = evalin('base','Process');
SeqControl = evalin('base','SeqControl');
SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
currentDepth = evalin('base','Receive(P.numRays+1).endDepth');
tof = ceil(2*currentDepth/demodFreq);
if SeqControl(6).argument < tof
    SeqControl(6).argument = tof;
    P.dopPRF = round(1/(tof*1e-06));
    SeqControl(6).argument = round(1/(P.dopPRF*1e-06));
    fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
    fprintf(['"dopPRF" is adjusted to ' num2str(P.dopPRF) '\n']);
    UI = evalin('base','UI');
    set(UI(7).handle(2),'Value',P.dopPRF);
    set(UI(7).handle(3),'String',num2str(P.dopPRF));
    for k = 1:2:length(Process(2).Parameters)
        if strcmp(Process(2).Parameters{k},'prf'), Process(2).Parameters{k+1} = P.dopPRF; end
    end
    assignin('base','Process',Process);
end
%---------------check Doppler PRF--------------------

assignin('base','P',P);
assignin('base','SeqControl',SeqControl);

% TGC
TGC = evalin('base','TGC');
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));
TGC(2).rangeMax = P.dopEndDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));
assignin('base','TGC',TGC);
assignin('base','nTGC',get(findobj('tag','TGCnum'),'Value'));
evalin('base','if VDAS==1, Result = loadTgcWaveform(nTGC); end');

Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',2,'prf',P.dopPRF};
Control(2).Command = 'update&Run';
Control(2).Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','SeqControl','TX','Receive','Recon','TGC'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');

return
%RangeChangeCallback


