%% This script will image the femoral artery in the pig
%% Plane Wave Bmode and Plane Wave Power Doppler
clear all
counter = 1;
%% BMode User Parameters
P.startDepthMm = 0;   % start depth in mm
P.endDepthMm = 25;    % end depth in mm
P.fovWidthMm = 15;    % the width of the bottom of the image (mm)
P.numAvgs = 1;        % number of averages for the Bmode image
P.numPlanes = 50;    % number of planes used to create a Bmode image
P.maxFPS = 10;        % Frame per second for Bmode imaging
P.c = 1540;           % Speed of Sound (m/s)      
P.beta = 4;         % Narrowness of the kaiser window
P.useElem = [1, 64];

%% Doppler User Parameters
dop.PRI = 68;        % Pulse repetition interval (number of images in an ensemble)
dop.maxV = 15; %cm/s   % Maximum velocity you're interested in measuring
dop.startDepthMm = 0; % mm
dop.endDepthMm = 15;  % mm
dop.vesDepth = 15;
dop.imWidthMm = 9;
dop.cycles = 4;
dop.numAngs = 10; %should be odd to allow for tx at 0 degrees
dop.vesR = 0.25; %mm
dop.prcPwr = 0;
dop.prcRF = 0;
dop.lastInterFrame = 0;
dop.outDir = 'FilterComparison';
dop.saveData = 0;
dop.fname = 'StringPhantom_BPFilt';
dop.cpl = 70;   
dop.vesAng = 45; %degrees

%% Specify system parameters.
Resource.Parameters.numTransmit = 64;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % number of receive channels.
Resource.Parameters.speedOfSound = P.c;
Resource.Parameters.verbose = 1;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.waitForProcessing = 0; % DMA transfers wait for processing.
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = 1;

%% Specify Trans structure array.
Trans.name = 'MUSIC_Biobox64';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans.frequency = 10.41666;
TWFreq = 250/floor(250/(4*Trans.frequency))/4;
Trans.frequency = TWFreq;
Trans.wvl = Resource.Parameters.speedOfSound*1000/(Trans.frequency*10^6);
Trans.Bandwidth = [7 14];
Trans = computeTrans2(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 50;  % set maximum high voltage limit for pulser supply.
Trans.name = 'custom';

dop.PRF = 4*dop.maxV/(P.c*100)*Trans.frequency*10^6;

%%
P.startDepth = P.startDepthMm/Trans.wvl;
P.endDepth = P.endDepthMm/Trans.wvl;
P.fovWidth = P.fovWidthMm/Trans.wvl;
ang = atan(-P.fovWidth/2/(P.endDepth));
P.angRange = linspace(-ang,ang,P.numPlanes);

dop.startDepth = dop.startDepthMm/Trans.wvl;
dop.endDepth = dop.endDepthMm/Trans.wvl;
dop.imWidth = dop.imWidthMm/Trans.wvl;

%% Specify PData structure array.
PData(1).PDelta = [Trans.spacing/2, 0, 1];  % x, y, z pdeltas
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil(P.fovWidth/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
PData(1).Origin = [-P.fovWidth/2,Trans.ElementPos(P.useElem(1),2),P.startDepth]; % x,y,z of upper lft crnr
PData(1).Coord = 'rectangular';

PData(1).Region = repmat(struct('Shape',struct('Name','Parallelogram',...
    'Position',[0,PData(1).Origin(2),0],...
    'width',Trans.elemInArr*Trans.spacing,...
    'height',P.endDepth-P.startDepth,...
    'angle',P.angRange(1))),1,P.numPlanes);

for i = 1:P.numPlanes
    PData(1).Region(i).Shape.angle = P.angRange(i);  
end
PData(1).Region = computeRegions(PData(1));

%% Power Doppler PData
maxAngDop = pi/2-atan((dop.vesDepth+dop.vesR)/(Trans.spacingMm*32));
if dop.numAngs ~= 1
    dop.dopAngs = linspace(-maxAngDop,maxAngDop,dop.numAngs);
else
    dop.dopAngs = 0;
end

PData(2).Coord = 'rectangular';
PData(2).PDelta = [Trans.spacing,0,dop.cycles];
PData(2).Size(1) = ceil((dop.endDepth-dop.startDepth)/PData(2).PDelta(3));
PData(2).Size(2) = ceil(dop.imWidth/PData(2).PDelta(1));
PData(2).Size(3) = 1;
PData(2).Origin = [-dop.imWidth/2, PData(1).Origin(2), dop.startDepth];
PData(2).Region = computeRegions(PData(2));

%% Specify Resources.
assert(1/(dop.PRF*dop.numAngs)>2*sqrt((Trans.elemInArr*Trans.spacing/2)^2+dop.endDepth^2)*Trans.wvl/(Resource.Parameters.speedOfSound*1000),'Plane Wave Imaging depth is too large for the given maxV');

maxAcqSamplesB = 128*ceil(2*sqrt((Trans.elemInArr*Trans.spacing)^2+P.endDepth^2)*5/128);
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = maxAcqSamplesB*P.numPlanes*1.5;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 2;     % 10 frames for RF data.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 2;

%Power Doppler and Background Ultrasound
dopInterBufFrames = 1;
maxAcqSamples = 128*ceil(2*sqrt((Trans.elemInArr*Trans.spacing)^2+dop.endDepth^2)*4/128);
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = maxAcqSamples*dop.numAngs*dop.PRI+maxAcqSamplesB*P.numPlanes*P.numAvgs;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 2;     % 10 frames for RF data.
Resource.InterBuffer(2).numFrames = 8;  % one intermediate buffer needed.
Resource.InterBuffer(2).pagesPerFrame = dop.PRI; 
Resource.InterBuffer(2).rowsPerFrame = PData(2).Size(1);
Resource.InterBuffer(2).colsPerFrame = PData(2).Size(2);
Resource.ImageBuffer(2).numFrames = 8;
Resource.ImageBuffer(2).rowsPerFrame = PData(2).Size(1);
Resource.ImageBuffer(2).colsPerFrame = PData(2).Size(2);

%Inter and Image Buffer is for processing background Bmode for Pwr Doppler
Resource.InterBuffer(3).numFrames = 1;  % one intermediate buffer needed.
Resource.InterBuffer(3).rowsPerFrame = PData(1).Size(1);
Resource.InterBuffer(3).colsPerFrame = PData(1).Size(2);
Resource.ImageBuffer(3).numFrames = 2;
Resource.ImageBuffer(3).rowsPerFrame = PData(1).Size(1);
Resource.ImageBuffer(3).colsPerFrame = PData(1).Size(2);

Resource.DisplayWindow(1).Title = 'MUSIC Bmode';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
%Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%               DwWidth, DwHeight];
Resource.DisplayWindow(1).Position = [250-DwWidth,(ScrnSize(4)-(DwHeight+550))/2, ...  
                    DwWidth, DwHeight];
                       
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(2),PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 400;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

Resource.DisplayWindow(2).Title = 'MUSIC Doppler';
Resource.DisplayWindow(2).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);
%Resource.DisplayWindow(2).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
%                                     DwWidth, DwHeight];
Resource.DisplayWindow(2).Position = [250+DwWidth*0.2,(ScrnSize(4)-(DwHeight+550))/2, ...  
                    DwWidth, DwHeight];

Resource.DisplayWindow(2).ReferencePt = [PData(1).Origin(1),PData(1).Origin(2),PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).numFrames = 5;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = grayscaleCFImap;
Resource.DisplayWindow(2).splitPalette = 1;

%% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

%Doppler Transmit
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,.67,dop.cycles*2,1];
%% Specify m TX structure arrays. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements),...
                   'TXPD',[],...
                   'peakBLMax',2,...
                   'peakCutOff',0.8,...
                   'aperture',1),1,(P.numPlanes+dop.numAngs));

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency*1000/(Resource.Parameters.speedOfSound);
end
for i = 1:P.numPlanes
    TX(i).Apod(P.useElem(1):P.useElem(2)) = 1;% kaiser(Trans.elemInArr,P.beta);
    TX(i).Origin = [0,Trans.ElementPos(P.useElem(1),2),0];
    TX(i).Steer = [P.angRange(i),0];
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(1));
    TX(i).peakBLMax = 32;
    TX(i).peakCutOff = 0.8;
end

lastBmodeTxArr1 = P.numPlanes;

for i = 1:dop.numAngs
    TX(lastBmodeTxArr1+i).waveform = 2;
    TX(lastBmodeTxArr1+i).Origin = [0,Trans.ElementPos(P.useElem(1),2),0];
    TX(lastBmodeTxArr1+i).Steer = [dop.dopAngs(i),0];
    TX(lastBmodeTxArr1+i).Apod(P.useElem(1):P.useElem(2)) = 1;%kaiser(Trans.elemInArr,P.beta);
%     TX(lastBmodeTxArr1+i).Steer = [0,0];
%     TX(lastBmodeTxArr1+i).focus = dop.focalDepth;
    TX(lastBmodeTxArr1+i).Delay = computeTXDelays(TX(lastBmodeTxArr1+i));
    TX(lastBmodeTxArr1+i).peakBLMax = 7;
    TX(lastBmodeTxArr1+i).peakCutOff = 0.5;
%     TX(lastBmodeTxArr1+i).TXPD = computeTXPD(TX(lastBmodeTxArr1+i),PData(2));
end

lastDopTxArr1 = lastBmodeTxArr1 + dop.numAngs;

%% TPC
TPC(1).name = 'US Power';
TPC(1).maxHighVoltage = 25;
TPC(1).hv = 1.6;

TPC(2).name = 'Doppler Power';
TPC(2).maxHighVoltage = 60;
TPC(2).hv = 5;
%%
% Specify TGC Waveform structure.
% TGC(1).CntrlPts = [189,314,457,698,770,911,948,976];
TGC(1).CntrlPts = ones(1,8)*1023;
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

TGC(2).CntrlPts = ones(1,8)*1023;
TGC(2).rangeMax = dop.endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));
%% Specify Receive structure arrays.
% - We need m Receive structures for each frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + P.fovWidth^2));
if abs(Trans.frequency-10.4167)<0.1
    coefbp = load('BP9_11MHz_fs=41.66668MHzv2.mat');
    bpfilt = coefbp.BP9_11MHz(1:ceil(length(coefbp.BP9_11MHz)/2));
%     bpfilt(1:2:end)=0;
    
    coeflp = load('LP14MHz_fs=41.666668MHzv2.mat');
    lpfilt = coeflp.LP14MHz(1:ceil(length(coeflp.LP14MHz)/2));
%     lpfilt(1:2:end) = 0;
%     bpfilt = [];
%     lpfilt = [];
else
    bpfilt = [];
    lpfilt = [];
end

Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','NS200BW', ...
                        'mode', 1, ...
                        'LowPassCoef', [],...
                        'InputFilter', [],...
                        'callMediaFunc', 0,...
                        'aperture',1), 1, Resource.RcvBuffer(1).numFrames*P.numAvgs*P.numPlanes + Resource.RcvBuffer(2).numFrames*(dop.numAngs*dop.PRI+P.numPlanes*P.numAvgs));%(Resource.RcvBuffer(1).numFrames*P.numAvgs*P.numPlanes + Resource.RcvBuffer(2).numFrames*(dop.numAngs*dop.PRI+P.numPlanes*P.numAvgs)+Resource.RcvBuffer(3).numFrames*(cdi.numLines*dop.PRI+P.numPlanes*P.numAvgs))
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numPlanes
        for k = 1:P.numAvgs
            if k == 1
                Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).mode = 0;
            end
            Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).Apod(P.useElem(1):P.useElem(2)) = 1;
            Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).framenum = i;
            Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).acqNum = (j-1)*P.numAvgs+k;
            Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).sampleMode = 'BS100BW';
        end
    end
end
lastBmodeRcvArr1 = Resource.RcvBuffer(1).numFrames*P.numPlanes*P.numAvgs;

for i = 1:Resource.RcvBuffer(2).numFrames
    currframe_rcv = lastBmodeRcvArr1+(i-1)*(P.numPlanes*P.numAvgs+dop.numAngs*dop.PRI);
    for k = 1:P.numPlanes
        for l = 1:P.numAvgs
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).bufnum = 2;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).mode = 0;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).Apod(P.useElem(1):P.useElem(2)) = 1;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).framenum = i;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).acqNum = (k-1)*P.numAvgs+l;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).endDepth = P.endDepth;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).startDepth = P.startDepth;
            Receive(currframe_rcv+(k-1)*P.numAvgs+l).sampleMode = 'BS100BW';

%             Receive(currframe_rcv+(k-1)*P.numAvgs+l).sampleMode = 'NS200BW';
        end
    end
    bkgdUSRcv = currframe_rcv+P.numPlanes*P.numAvgs;
    
    for j = 1:dop.numAngs*dop.PRI
        Receive(bkgdUSRcv+j).bufnum = 2;
        Receive(bkgdUSRcv+j).mode = 0;
        currBeam = mod(j,dop.numAngs);
        Receive(bkgdUSRcv+j).Apod(P.useElem(1):P.useElem(2)) = 1; 
        Receive(bkgdUSRcv+j).framenum = i;
        Receive(bkgdUSRcv+j).acqNum = P.numPlanes*P.numAvgs+j;
        Receive(bkgdUSRcv+j).endDepth = dop.endDepth;
        Receive(bkgdUSRcv+j).startDepth = dop.startDepth;
        Receive(bkgdUSRcv+j).sampleMode = 'NS200BW';
        Receive(bkgdUSRcv+j).LowPassCoef = lpfilt;
        Receive(bkgdUSRcv+j).InputFilter = bpfilt;

    end
end
lastDopRcvArr1 = lastBmodeRcvArr1 + Resource.RcvBuffer(2).numFrames*(dop.numAngs*dop.PRI+P.numPlanes*P.numAvgs);

%% RcvProfile
RcvProfile(1).AntiAliasCutoff = 20;
RcvProfile(1).PgaGain = 30;
RcvProfile(1).LnaGain = 24;
RcvProfile(1).LnaHPF = 200;
RcvProfile(1).LnaZinSel = 28;

RcvProfile(2).AntiAliasCutoff = 20;
RcvProfile(2).PgaGain = 30; %30
RcvProfile(2).LnaGain = 24;
RcvProfile(2).LnaHPF = 200;
RcvProfile(2).LnaZinSel = 28;

%% Specify Recon structure arrays.
% - We need one Recon structure.  Each frame will use
%   m ReconInfo structures, since we are using m
%   synthetic aperture acquisitions.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', []),1,3);

% Recon(2).newFrameTimeout = 1000;

           
%Doppler
Recon(2).pdatanum = 2;
Recon(2).rcvBufFrame = -1;
Recon(2).IntBufDest = [2,-1];
Recon(2).ImgBufDest = [0, 0];
Recon(2).senscutoff = 0.7;
% Recon(2).newFrameTimeout = 10000;

%Background Bmode for Pwr Doppler
Recon(3).pdatanum = 1;
Recon(3).rcvBufFrame = -1;
Recon(3).IntBufDest = [3, 1];
Recon(3).ImgBufDest = [3, -1];
Recon(3).senscutoff = 0.7;
%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, (P.numPlanes+dop.numAngs*dop.PRI+P.numPlanes));%
% - Set specific ReconInfo attributes.

ReconInfo(1).Pre = 'clearInterBuf';
ReconInfo(1).mode = 'replaceIQ';  % replace IQ data
for j = 1:P.numPlanes 
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
% ReconInfo(P.numPlanes).mode = 'accumIQ_replaceIntensity'; % accum & detect
ReconInfo(P.numPlanes).Post = 'IQ2IntensityImageBuf'; % accum & detect
lastBmodeRIArr1 = P.numPlanes;
Recon(1).RINums = 1:lastBmodeRIArr1;

%Pwr Doppler Arr 1 RI
ReconInfo(lastBmodeRIArr1+1).Pre = 'clearInterBuf';
ReconInfo(lastBmodeRIArr1+1).mode = 'replaceIQ';
for k = 1:dop.PRI
    ReconInfo(lastBmodeRIArr1+(k-1)*dop.numAngs+1).mode = 'replaceIQ';
    for j = 1:dop.numAngs  % For each row in the column
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numAngs+j).txnum = lastBmodeTxArr1+j;
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numAngs+j).rcvnum = lastBmodeRcvArr1+P.numPlanes*P.numAvgs+(k-1)*dop.numAngs+j;
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numAngs+j).pagenum = k;
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numAngs+j).regionnum = 1;
    end
end
lastDopRIArr1 = lastBmodeRIArr1 + dop.numAngs*dop.PRI;
Recon(2).RINums = lastBmodeRIArr1+1:lastDopRIArr1;

%Background Ultrasound for Pwr Doppler image
ReconInfo(lastDopRIArr1+1).Pre = 'clearInterBuf';
ReconInfo(lastDopRIArr1+1).mode = 'replaceIQ';  % replace IQ data
for j = 1:P.numPlanes
    for k = 1:P.numAvgs
        ReconInfo(lastDopRIArr1+j).txnum = j;
        ReconInfo(lastDopRIArr1+j).rcvnum = lastBmodeRcvArr1+(j-1)*P.numAvgs+k;
        ReconInfo(lastDopRIArr1+j).regionnum = j;
    end
end
ReconInfo(lastDopRIArr1+P.numPlanes).Post = 'IQ2IntensityImageBuf'; % accum & detect
lastDopBmodeRIArr1 = lastDopRIArr1+P.numPlanes;
Recon(3).RINums = lastDopRIArr1+1:lastDopBmodeRIArr1;
%% Specify Process structure array.
pers = 20;
dopThresh = 0.01;
maxPwr = 20;
compFactorB = 40;
rejB = 0;
pGainB = 2;
dispThresh = 1;

compFactorD = 20;
rejD = 0;
pGainD = 3;
bmodeProcess = 1;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',pGainB,...            % pgain is image processing gain
                         'reject',rejB,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','reduceSpeckle2',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',compFactorB,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
PDIProcess = 2;
Process(2).classname = 'Doppler';
Process(2).method = 'computeCFIPowerEst';
Process(2).Parameters = {'IntBufSrc',[2,-1],...   
                         'SrcPages',[4 dop.PRI-8],...   
                         'ImgBufDest',[2,-1],...            
                         'pdatanum',2,...
                         'prf',ceil(dop.PRF),...
                         'wallFilter','WeakFlowLow',...
                         'pwrThreshold',dopThresh,...
                         'maxPower',maxPwr};
                     
%Display pwr Doppler Arr 1
dopDispProcess = 3;
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum', -1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',pGainD,...            % pgain is image processing gain
                         'reject',rejD,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','power',...
                         'compressFactor',compFactorD,...
                         'mappingMethod','upperHalf',...
                         'display',1,...      % display image after processing
                         'displayWindow',2,...
                         'srcData','unsignedColor',...
                         'threshold',dop.cpl};

bkgdBmodeProcess = 4;
%Display Bkgd Ultrasound
Process(4).classname = 'Image';
Process(4).method = 'imageDisplay';
Process(4).Parameters = {'imgbufnum',3,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',pGainB,...            % pgain is image processing gain
                         'reject',rejB,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','reduceSpeckle2',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',compFactorB,...
                         'mappingMethod','lowerHalf',...
                         'display',0,...      % display image after processing
                         'displayWindow',2};

% %Display Bkgd Ultrasound
% Process(8).classname = 'Image';
% Process(8).method = 'imageDisplay';
% Process(8).Parameters = {'imgbufnum',3,...   % number of buffer to process.
%                          'framenum',-1,...   % (-1 => lastFrame)
%                          'pdatanum',1,...    % number of PData structure to use
%                          'pgain',pGainB,...            % pgain is image processing gain
%                          'reject',rejB,...      % reject level
%                          'persistMethod','simple',...
%                          'persistLevel',pers,...
%                          'interpMethod','4pt',...
%                          'grainRemoval','none',...
%                          'processMethod','reduceSpeckle2',...
%                          'averageMethod','none',...
%                          'compressMethod','power',...
%                          'compressFactor',compFactorB,...
%                          'mappingMethod','lowerHalf',...
%                          'display',0,...      % display image after processing
%                          'displayWindow',2};
% 
% %Display Doppler Arr 2
% Process(9).classname = 'Image';
% Process(9).method = 'imageDisplay';
% Process(9).Parameters = {'imgbufnum',4,...   % number of buffer to process.
%                          'framenum',-1,...   % (-1 => lastFrame)
%                          'pdatanum',2,...    % number of PData structure to use
%                          'pgain',pGainD,...            % pgain is image processing gain
%                          'reject',rejD,...      % reject level
%                          'persistMethod','simple',...
%                          'persistLevel',pers,...
%                          'interpMethod','4pt',...
%                          'grainRemoval','none',...
%                          'processMethod','reduceSpeckle2',...
%                          'averageMethod','none',...
%                          'compressMethod','power',...
%                          'compressFactor',compFactorD,...
%                          'mappingMethod','lowerHalf',...
%                          'display',1,...      % display image after processing
%                          'displayWindow',2};
%                      
% 
% Process(5).classname = 'External';
% Process(5).method = 'myTic';
% Process(5).Parameters = {'srcbuffer','none',...
%                          'dstbuffer','none'};
%                      
%                      

%                      
% %Display CDI

%                      
% Process(9).classname = 'External';
% Process(9).method = 'calcMeanVel';
% Process(9).Parameters = {'srcbuffer','inter',...
%                          'srcbufnum', 2,...
%                          'srcframenum', -1,...
%                          'dstbuffer','image',...
%                          'dstbufnum',2,...
%                          'dstframenum',-1};
% 
% Process(10).classname = 'External';
% Process(10).method = 'dopplerRecon';
% Process(10).Parameters = {'srcbuffer','receive',...
%                          'srcbufnum', 2,...
%                          'srcframenum', -1,...
%                          'dstbuffer','inter',...
%                          'dstbufnum',2,...
%                          'dstframenum',-1};
%                      
% Process(11).classname = 'External';
% Process(11).method = 'saveIQData';
% Process(11).Parameters = {'srcbuffer','inter',...
%                          'srcbufnum', 2,...
%                          'srcframenum', -1,...
%                          'dstbuffer','none'};

%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
% SeqControl(1).argument = 1;

SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
maxdelay = ceil(2*sqrt((Trans.ElementPos(64)*Trans.wvl*2)^2 + P.endDepthMm^2)/(Resource.Parameters.speedOfSound*1000)*1E6*1.5);
SeqControl(2).argument = ceil(1/(P.maxFPS*P.numPlanes)*1E6);  % 200 usec

SeqControl(3).command = 'returnToMatlab';

SeqControl(4).command = 'sync';
SeqControl(4).argument = 5000000;

SeqControl(5).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(5).argument = 5000;  % 200 usec

sc_dop_ttna = 6;
SeqControl(6).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(6).argument = ceil(1/(dop.PRF*dop.numAngs)*1E6);

SeqControl(7).command = 'setRcvProfile';  % time between synthetic aperture acquisitions
SeqControl(7).argument = 1; 

SeqControl(8).command = 'setTPCProfile';  % time between synthetic aperture acquisitions
SeqControl(8).argument = 1;
SeqControl(8).condition = 'immediate';

SeqControl(9).command = 'noop';  % time between synthetic aperture acquisitions
SeqControl(9).argument = 80000/0.2; %8 ms

SeqControl(10).command = 'setRcvProfile';  % time between synthetic aperture acquisitions
SeqControl(10).argument = 2; 

SeqControl(11).command = 'jump'; % jump back to start.

SeqControl(12).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(12).argument = 100000;  % 200 usec

SeqControl(13).command = 'setTPCProfile';  % time between synthetic aperture acquisitions
SeqControl(13).argument = 2;
SeqControl(13).condition = 'immediate';

SeqControl(14).command = 'returnToMatlab'; % jump back to start.

SeqControl(15).command = 'triggerOut';

SeqControl(16).command = 'setTPCProfile';  % time between synthetic aperture acquisitions
SeqControl(16).argument = 3;
SeqControl(16).condition = 'immediate';

SeqControl(17).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(17).argument = 10;  % 200 usec

SeqControl(18).command = 'jump'; % jump back to start.
SeqControl(18).argument = 1;  % 200 usec

SeqControl(19).command = 'sync';
SeqControl(19).argument =100;

SeqControl(20).command = 'jump'; % jump back to start.
SeqControl(20).argument = 1;  % 200 usec

SeqControl(21).command = 'timeToNextAcq';
SeqControl(21).argument = ceil(1/(dop.PRF)*1E6);

SeqControl(22).command = 'timeToNextAcq';
SeqControl(22).argument = maxdelay;

nsc = 23; % nsc is count of SeqControl objects

%% Specify Event structure array.
n = 1;
Event(n).info = 'Set Rcv Profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 7;
n = n+1;

Event(n).info = 'Set TPC';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 8;
n = n+1;

Event(n).info = 'No op';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 9;
n = n+1;

startBMode = n;
SeqControl(11).argument = startBMode;
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numPlanes
        for k = 1:P.numAvgs
            Event(n).info = 'Plane Wave TX';
            Event(n).tx = j;
            Event(n).rcv = (i-1)*P.numAvgs*P.numPlanes+(j-1)*P.numAvgs+k;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 2;
            n = n+1;
        end
    end
    Event(n-1).seqControl = [5,nsc]; % modify last acquisition Event's seqControl
        SeqControl(nsc).command = 'transferToHost'; % transfer frame to host buffer
        nsc = nsc+1;
    
    Event(n).info = 'Recon and Display';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;
 
end

Event(n).info = 'Return to Matlab';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 3;
n = n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 11;
n = n+1;

startDopArr1 = n;
Event(n).info = 'Set Pwr Doppler Rcv Profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 10;
n = n+1;

Event(n).info = 'Set Pwr Doppler TPC Profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 13;
n = n+1;

Event(n).info = 'No op';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 9;
n = n+1;

SeqControl(1).argument = n;
for l = 1:Resource.RcvBuffer(2).numFrames
    currframe = (l-1)*(P.numPlanes*P.numAvgs+dop.numAngs*dop.PRI);
    for j = 1:P.numPlanes
        for k = 1:P.numAvgs
            Event(n).info = 'Bkgd Plane Wave TX';
            Event(n).tx = j;
            Event(n).rcv = lastBmodeRcvArr1+currframe+(j-1)*P.numAvgs+k;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 22;
            n = n+1;
        end
    end

    for i = 1:dop.PRI
        for j = 1:dop.numAngs
            Event(n).info = 'Pwr Doppler Plane TX';
            Event(n).tx = lastBmodeTxArr1 + j;
            Event(n).rcv = lastBmodeRcvArr1 + currframe + P.numPlanes*P.numAvgs + (i-1)*dop.numAngs + j;
            Event(n).recon = 0;
            Event(n).process = 0;
            Event(n).seqControl = 6;
            n = n+1;
        end
    end
    Event(n-1).seqControl = [nsc,12];
        SeqControl(nsc).command = 'transferToHost';
        nsc = nsc + 1;

    Event(n).info = 'Recon';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = [3,2];
    Event(n).process = 2;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Display Bkgd Bmode';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 4;
    Event(n).seqControl = 0;
    n = n+1;

    Event(n).info = 'Display Doppler';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;
end
Event(n).info = 'Return to MATLAB';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 14;
n = n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

%% UI
import vsv.seq.uicontrol.VsButtonGroupControl

UI(1).Control = {'UserB6','Style','VsToggleButton','Label','PDI'};
UI(1).Callback = text2cell('%SingleAcqDoppler');

UI(2).Control =  {'UserB4','Style','VsSlider','Label','PDI Priority',...
                  'SliderMinMaxVal',[0,255,dop.cpl],...
                  'SliderStep',[1/255,0.1],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%changePriority');

UI(3).Control =  {'UserB3','Style','VsSlider','Label','PRF (kHz)',...
                  'SliderMinMaxVal',[0,19,dop.PRF/1000],...
                  'SliderStep',[0.1,1],'ValueFormat','%0.1f'};
UI(3).Callback = text2cell('%changePRF');

UI(4).Control =  {'UserB2','Style','VsSlider','Label','BMode Reject',...
                  'SliderMinMaxVal',[0,100,rejB],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%2d'};
UI(4).Callback = text2cell('%bmodereject');

UI(5).Control =  {'UserC1','Style','VsSlider','Label','BMode Compress',...
                  'SliderMinMaxVal',[0,1,compFactorB/100],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%0.2f'};
UI(5).Callback = text2cell('%changeCompression');

UI(6).Control =  {'UserB1','Style','VsSlider','Label','BMode Gain',...
                  'SliderMinMaxVal',[0,100,pGainB],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%2d'};
UI(6).Callback = text2cell('%changeGain');

UI(7).Control =  {'UserC7','Style','VsSlider','Label','Disp Thresh',...
                  'SliderMinMaxVal',[1,101,dispThresh],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%2d'};
UI(7).Callback = text2cell('%changeDispThresh');

UI(8).Control =  {'UserB5','Style','VsSlider','Label','PDI Pwr Threshold',...
                  'SliderMinMaxVal',[0,1,dopThresh],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%0.2f'};
UI(8).Callback = text2cell('%changePwrThresh');

UI(9).Control = {'UserB7','Style','VsButtonGroup','Title','Doppler Mode','NumButtons',2,'Labels',{'Velocity','Power'}};
UI(9).Callback = text2cell('%DopplerMode');

% Save all the structures to a .mat file.
filename = 'C:\Users\Administrator\Documents\Vantage-4.8.4-2211151000\MatFiles\MUSIC_PlaneDoppler.mat';
save(filename);
VSX;

if dop.saveData == 1
    t = yyyymmdd(datetime);
    if ~exist(dop.outDir,'dir')
        mkdir(dop.outDir);
    end
    fname = [dop.outDir,'\',sprintf('%d_%s.mat',t,dop.fname)];
    saveDopData(fname)
end

return

% **** Callback routines to be converted by text2cell function. ****
%% UI(1) - Single Acquisition Doppler
%SingleAcqDoppler
singleshot_start = evalin('base','startDopArr1');
Control = evalin('base','Control');
Control.Command = 'set&Run';
if UIState
    Control.Parameters = {'Parameters',1,'startEvent',singleshot_start};
else
    Control.Parameters = {'Parameters',1,'startEvent',1};
end
assignin('base','Control', Control);

%SingleAcqDoppler

%% UI(4) - change color priority
%changePriority

dop = evalin('base','dop');
dopDispProcess = evalin('base','dopDispProcess');

Process = evalin('base','Process');
for k = 1:2:length(Process(dopDispProcess).Parameters)
    if strcmp(Process(dopDispProcess).Parameters{k},'threshold'), Process(dopDispProcess).Parameters{k+1} = UIValue; end
end

Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',dopDispProcess,'threshold',UIValue};
dop.cpl = UIValue;
assignin('base','dop',dop);
assignin('base','Control', Control);

%changePriority

%% UI(5) - Change PRF
%changePRF
dop = evalin('base','dop');
PDIProcess = evalin('base','PDIProcess');
sc_dop_ttna = evalin('base','sc_dop_ttna');
SeqControl = evalin('base','SeqControl');
Resource = evalin('base','Resource');
startBMode = evalin('base','startBMode');
lastBmodeRcvArr1 = evalin('base','lastBmodeRcvArr1');

Receive = evalin('base','Receive');
fs = Receive(lastBmodeRcvArr1).decimSampleRate/Receive(lastBmodeRcvArr1).quadDecim;
dop.PRF = UIValue*1000;
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',PDIProcess,'prf',dop.PRF};

setLoopCntTo1 = evalin('base','setLoopCntTo1');
SeqControl(sc_dop_ttna).argument = 1/(dop.numAngs*dop.PRF*1E-6);
assignin('base','SeqControl',SeqControl);
Control(2).Command = 'update&Run';
Control(2).Parameters = {'SeqControl'};
Control(3).Command = 'set&Run';
Control(3).Parameters = {'Parameters',1,'startEvent',setLoopCntTo1};

disp([num2str(dop.PRI*((dop.numAngs)*SeqControl(sc_dop_ttna).argument)/1000),...
    ' ms singleshot acquisition time; ', num2str(dop.PRF), ' kHz effective PRF; ',...
    num2str(100*(Resource.Parameters.speedOfSound*dop.PRF)/(fs*10^6)),...
    ' cm/sec max sampled velocity'])
dop.maxV = 100*(Resource.Parameters.speedOfSound*dop.PRF)/(fs*10^6)/4; % cm/sec
assignin('base','dop',dop);
assignin('base','Control', Control);

%changePRF

%bmodereject
rejB = evalin('base','rejB');
bkgdBmodeProcess = evalin('base','bkgdBmodeProcess');
bmodeProcess = evalin('base','bmodeProcess');
rejB = UIValue;
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',bmodeProcess,'reject',rejB};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',bkgdBmodeProcess,'reject',rejB};
assignin('base','Control',Control);
assignin('base','rejB',rejB);
%bmodereject

%changeCompression
compFactorB = evalin('base','compFactorB');
bkgdBmodeProcess = evalin('base','bkgdBmodeProcess');
bmodeProcess = evalin('base','bmodeProcess');
compFactorB = UIValue;
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',bmodeProcess,'compression',compFactorB};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',bkgdBmodeProcess,'compression',compFactorB};
assignin('base','Control',Control);
assignin('base','compFactorB',compFactorB);
%changeCompression

%changeGain
pGainB = evalin('base','pGainB');
bkgdBmodeProcess = evalin('base','bkgdBmodeProcess');
bmodeProcess = evalin('base','bmodeProcess');
pGainB = UIValue;
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',bmodeProcess,'pgain',pGainB};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',bkgdBmodeProcess,'pgain',pGainB};
assignin('base','Control',Control);
assignin('base','pGainB',pGainB);
%changeGain

%changeDispThresh
dispThresh = evalin('base','dispThresh');
bkgdBmodeProcess = evalin('base','bkgdBmodeProcess');
Process = evalin('base','Process');
Control = evalin('base','Control');

for k = 1:2:length(Process(bkgdBmodeProcess).Parameters)
    if strcmp(Process(bkgdBmodeProcess).Parameters{k},'threshold'), Process(bkgdBmodeProcess).Parameters{k+1} = UIValue; end
end

Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',bkgdBmodeProcess,'threshold',UIValue};

dispThresh = UIValue;
assignin('base','Process',Process);

assignin('base','Control',Control);
assignin('base','dispThresh',dispThresh);

%changeDispThresh

%% UI(3) - change power threshold
%changePwrThresh
dop = evalin('base','dop');
PDIProcess = evalin('base','PDIProcess');

Process = evalin('base','Process');
for k = 1:2:length(Process(PDIProcess).Parameters)
    if strcmp(Process(PDIProcess).Parameters{k},'pwrThreshold'), Process(PDIProcess).Parameters{k+1} = UIValue; end
end

assignin('base','Process',Process);

% Set Control.Command to set Doppler threshold.
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',PDIProcess,'pwrThreshold',UIValue};

pwrThresh = UIValue;
assignin('base','pwrThresh',pwrThresh);
assignin('base','Control', Control);

%changePwrThresh

%DopplerMode
Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
Process = evalin('base','Process');
Resource = evalin('base','Resource');
dopProc = evalin('base','PDIProcess');
dopDisp = evalin('base','dopDispProcess');
pers = evalin('base','pers');

switch UIState
    
   case 1  % Velocity mode
      newMap = grayscaleCFImap;
      newMap(1:128,:) = Resource.DisplayWindow(2).Colormap(1:128,:);
      Resource.DisplayWindow(2).Colormap = newMap;
%       assignin('base','persp',get(findobj('Tag','UserB1Slider'),'Value'));

      Control(1).Parameters = {'Process',dopProc,'method','computeCFIFreqEst'};
      Control(2).Parameters = {'Process',dopDisp,'srcData','signedColor','persistMethod','dynamic','persistLevel',pers};
      Control(3).Parameters = {'DisplayWindow',2,'colormap',newMap};
      Control(4).Parameters = {'ImageBuffer',2,'lastFrame',0};
%       set(findobj('tag','UserB1Edit'),'String',num2str(persValue,'%3.0f'));
%       set(findobj('tag','UserB1Slider'),'Value',persValue);
      assignin('base','DopState','freq');
      % Set modified Process attributes in base Matlab environment.
      Process(dopProc).method = 'computeCFIFreqEst';
      for k = 1:2:length(Process(dopProc).Parameters)
          if strcmp(Process(dopProc).Parameters{k},'srcData'), Process(dopProc).Parameters{k+1} = 'signedColor';
          elseif strcmp(Process(dopProc).Parameters{k},'persistMethod'), Process(dopProc).Parameters{k+1} = 'dynamic';
          elseif strcmp(Process(dopProc).Parameters{k},'persistLevel'), Process(dopProc).Parameters{k+1} = pers;
          end
      end
   case 2  % Power mode
      newMap = grayscaleCPAmap;
      newMap(1:128,:) = Resource.DisplayWindow(2).Colormap(1:128,:);
      Resource.DisplayWindow(2).Colormap = newMap;
      for k = 1:2:length(Process(dopDisp).Parameters)
          if strcmp(Process(dopDisp).Parameters{k},'persistLevel'), pers = Process(dopDisp).Parameters{k+1}; end
      end
%       assignin('base','persf',persf);
%       persValue = evalin('base','persp');
      Control(1).Parameters = {'Process',dopProc,'method','computeCFIPowerEst'};
      Control(2).Parameters = {'Process',dopDisp,'srcData','unsignedColor','persistMethod','simple','persistLevel',pers};
      Control(3).Parameters = {'DisplayWindow',2,'colormap',newMap};
      Control(4).Parameters = {'ImageBuffer',2,'lastFrame',0};
%       set(findobj('tag','UserB1Edit'),'String',num2str(pers,'%3.0f'));
%       set(findobj('tag','UserB1Slider'),'Value',persValue);
      assignin('base','DopState','power');
      Process(dopProc).method = 'computeCFIPowerEst';
      for k = 1:2:length(Process(3).Parameters)
          if strcmp(Process(dopDisp).Parameters{k},'srcData'), Process(dopDisp).Parameters{k+1} = 'unsignedColor';
          elseif strcmp(Process(dopDisp).Parameters{k},'persistMethod'), Process(dopDisp).Parameters{k+1} = 'simple';
          elseif strcmp(Process(dopDisp).Parameters{k},'persistLevel'), Process(dopDisp).Parameters{k+1} = pers;
          end
      end
end

assignin('base','newMap',newMap);
evalin('base','Resource.DisplayWindow(1).Colormap = newMap;');
assignin('base','Process',Process);
assignin('base','Control', Control);

% If PTool window is open, adjust all uicontrols
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    posPTool = get(hPTool,'position');
    PTool;
    set(findobj('tag','ProcessTool'),'position',posPTool);
end

% If ColorMapTool is open, close it.
hCMTool = findobj('tag','ColorMapTool');
if ishandle(hCMTool),
    delete(hCMTool);
    set(findobj('tag','toolsMenu'),'Value',1); % set tools selection back to none
end

return
%DopplerMode