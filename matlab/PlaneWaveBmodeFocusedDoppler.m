%% HALEY USE THIS SCRIPT
%% Simple Plane Wave Transmit Script
%This script will fire from a single element and visualize the resulting
%RF data

%This is designed to be used with only 1 connector

% Start with Orange Board A connected to Green Connector L
% Collect Doppler Data at one position on the spinal cord
% Close GUI and type in the filename and the board that was connected to L
% Take screenshots of the doppler and US images and save them in folder 20230228_PigSurgery
% Repeat process with Board A connected to Connector H

clear all
useArr1 = 1;
%% BMode User Parameters
P.startDepthMm = 0;   % start depth in mm
P.endDepthMm = 25;    % end depth in mm
P.fovWidthMm = 15;    % the width of the bottom of the image (mm)
P.numAvgs = 1;        % number of averages for the Bmode image
P.numPlanes = 100;    % number of planes used to create a Bmode image
P.maxFPS = 50;        % Frame per second for Bmode imaging
P.c = 1540;           % Speed of Sound (m/s)      
P.beta = 5.6;         % Narrowness of the kaiser window
%% Doppler User Parameters
dop.PRI = 200;        % Pulse repetition interval (number of images in an ensemble)
dop.maxV = 5; %cm/s   % Maximum velocity you're interested in measuring
dop.startDepthMm = 5; % mm
dop.endDepthMm = 20;  % mm
dop.imWidthMm = 5;
dop.focalDepthMm = 11;      % Number of Angles used for Doppler imaging
dop.cycles = 8;
dop.fNum = 2;
dop.numBeams = 16;
%% Specify system parameters.
Resource.Parameters.numTransmit = 128;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 128;  % number of receive channels.
Resource.Parameters.speedOfSound = P.c;
Resource.Parameters.verbose = 1;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.waitForProcessing = 0; % DMA transfers wait for processing.
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = [1,2];

%% Specify Trans structure array.
Trans.name = 'MUSIC2';
Trans.units = 'mm'; % required in Gen3 to prevent default to mm units
Trans = computeTrans2(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 30;  % set maximum high voltage limit for pulser supply.
Trans.name = 'custom';
Trans.frequency = 10.4;
TWFreq = 250/floor(250/(4*Trans.frequency))/4;
Trans.frequency = TWFreq;
Trans.wvl = Resource.Parameters.speedOfSound*1000/(Trans.frequency*10^6);

dop.PRF = 4*dop.maxV/(P.c*100)*Trans.frequency*10^6;

%%
P.startDepth = P.startDepthMm/Trans.wvl;
P.endDepth = P.endDepthMm/Trans.wvl;
P.fovWidth = P.fovWidthMm/Trans.wvl;
% quarterArray = -Trans.spacing*Trans.elemInArr/4;
ang = atan(-P.fovWidth/2/(P.endDepth-P.startDepth));
angRange = linspace(-ang,ang,P.numPlanes);

dop.startDepth = dop.startDepthMm/Trans.wvl;
dop.endDepth = dop.endDepthMm/Trans.wvl;
dop.imWidth = dop.imWidthMm/Trans.wvl;
dop.focalDepth = dop.focalDepthMm/Trans.wvl;
%% Specify PData structure array.
PData(1).PDelta = [Trans.spacing, 0, 1];  % x, y, z pdeltas
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(2) = ceil(P.fovWidth/PData(1).PDelta(1));
PData(1).Size(3) = 1;      % single image page
if useArr1
    PData(1).Origin = [-P.fovWidth/2,0,P.startDepth]; % x,y,z of upper lft crnr
else
    PData(1).Origin = [-P.fovWidth/2,Trans.ElementPos(Trans.elemInArr+1,2),P.startDepth]; % x,y,z of upper lft crnr
end
PData(1).Region = repmat(struct('Shape',struct('Name','Parallelogram',...
    'Position',[0,PData(1).Origin(2),0],...
    'width',Trans.elemInArr*Trans.spacing,...
    'height',P.endDepth-P.startDepth,...
    'angle',angRange(1))),1,P.numPlanes);
for i = 1:P.numPlanes
    PData(1).Region(i).Shape.angle = angRange(i);  
end
PData(1).Region = computeRegions(PData(1));

%Doppler PData
centerFocus = linspace(-dop.imWidth/2,dop.imWidth/2,dop.numBeams);
elemInAp = ceil((dop.focalDepth/dop.fNum)/Trans.spacing);

if elemInAp > Trans.elemInArr
    elemInAp = Trans.elemInArr;
end

PData(2).PDelta = PData(1).PDelta;  % x, y, z pdeltas
PData(2).Size(1) = ceil((dop.endDepth-dop.startDepth)/PData(2).PDelta(3));
PData(2).Size(2) = ceil(dop.imWidth/PData(2).PDelta(1));
PData(2).Size(3) = 1;      % single image page
if useArr1
    PData(2).Origin = [-dop.imWidth/2,0,dop.startDepth]; % x,y,z of upper lft crnr
else
    PData(2).Origin = [-dop.imWidth/2,Trans.ElementPos(Trans.elemInArr+1,2)/Trans.wvl,dop.startDepth]; % x,y,z of upper lft crnr
end

PData(2).Region = repmat(struct('Shape',struct('Name','Parallelogram',...
    'Position',[0,PData(2).Origin(2),dop.startDepth],...
    'width',dop.imWidth,...
    'height',dop.endDepth-dop.startDepth,...
    'angle',0)),1,dop.numBeams);

for i = 1:dop.numBeams
    PData(2).Region(i).Shape.Position(1) = centerFocus(i);
    PData(2).Region(i).Shape.width = centerFocus(2)-centerFocus(1);
end

PData(2).Region = computeRegions(PData(2));

%% Specify Resources.
assert(1/(dop.PRF*dop.numBeams)>2*sqrt((Trans.elemInArr*Trans.spacing)^2+P.endDepth^2)*Trans.wvl/(Resource.Parameters.speedOfSound*1000),'Imaging depth is too large for the given maxV');

maxAcqSamples = 128*ceil(2*sqrt((Trans.elemInArr*Trans.spacing)^2+P.endDepth^2)*2.5/128);
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = maxAcqSamples*P.numPlanes;
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = 6;     % 10 frames for RF data.
Resource.InterBuffer(1).numFrames = 1;  % one intermediate buffer needed.
Resource.ImageBuffer(1).numFrames = 10;

Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = maxAcqSamples*dop.numBeams*dop.PRI;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = 6;     % 10 frames for RF data.
Resource.InterBuffer(2).numFrames = 6;  % one intermediate buffer needed.
Resource.InterBuffer(2).pagesPerFrame = dop.PRI; 
Resource.ImageBuffer(2).numFrames = 6;

% Resource.RcvBuffer(3).datatype = 'int16';
% Resource.RcvBuffer(3).rowsPerFrame = maxAcqSamples*P.numPlanes;
% Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
% Resource.RcvBuffer(3).numFrames = 6;     % 10 frames for RF data.
% Resource.InterBuffer(3).numFrames = 1;  % one intermediate buffer needed.
% Resource.ImageBuffer(3).numFrames = 10;
% 
% Resource.RcvBuffer(4).datatype = 'int16';
% Resource.RcvBuffer(4).rowsPerFrame = maxAcqSamples*dop.numAngs*dop.PRI;
% Resource.RcvBuffer(4).colsPerFrame = Resource.Parameters.numRcvChannels;
% Resource.RcvBuffer(4).numFrames = 6;     % 10 frames for RF data.
% Resource.InterBuffer(4).numFrames = 10;  % one intermediate buffer needed.
% Resource.ImageBuffer(4).numFrames = 10;
Resource.InterBuffer(3).datatype = 'complex single';
Resource.InterBuffer(3).numFrames = 10;  % one intermediate buffer needed.
Resource.ImageBuffer(3).numFrames = 4;


Resource.DisplayWindow(1).Title = 'MUSIC Bmode';
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),PData(1).Origin(2),PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Verasonics';
Resource.DisplayWindow(1).numFrames = 20;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

Resource.DisplayWindow(2).Title = 'MUSIC CDI';
Resource.DisplayWindow(2).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(2).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(2).ReferencePt = [PData(1).Origin(1),PData(1).Origin(2),PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(2).Type = 'Verasonics';
Resource.DisplayWindow(2).numFrames = 50;
Resource.DisplayWindow(2).AxesUnits = 'mm';
Resource.DisplayWindow(2).Colormap = grayscaleCFImap;
Resource.DisplayWindow(2).splitPalette = 1;

Resource.DisplayWindow(3).Title = 'MUSIC PDI';
Resource.DisplayWindow(3).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(2).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(2).pdelta);
Resource.DisplayWindow(3).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(3).ReferencePt = [PData(1).Origin(1),PData(1).Origin(2),PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(3).Type = 'Verasonics';
Resource.DisplayWindow(3).numFrames = 50;
Resource.DisplayWindow(3).AxesUnits = 'mm';
Resource.DisplayWindow(3).Colormap = grayscaleCPAmap;
Resource.DisplayWindow(3).splitPalette = 1;

%% Specify Transmit waveform structure.
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,.67,2,1];

%Doppler Transmit
TW(2).type = 'parametric';
TW(2).Parameters = [Trans.frequency,.67,dop.cycles,1];

%% Specify m TX structure arrays. Transmit on element n in the array for event n.
TX = repmat(struct('waveform', 1, ...
                   'Origin', zeros(1,3), ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Apod', zeros(1,Trans.numelements), ...
                   'Delay', zeros(1,Trans.numelements),...
                   'TXPD',[],...
                   'peakBLMax',2,...
                   'peakCutOff',2.3),1,(P.numPlanes+dop.numBeams));

scaleToWvl = 1;
if strcmp(Trans.units, 'mm')
    scaleToWvl = Trans.frequency*1000/(Resource.Parameters.speedOfSound);
end
for i = 1:P.numPlanes
    if useArr1 
        TX(i).Apod(1:Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
    else
        TX(i).Apod(Trans.elemInArr+1:2*Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
        TX(i).Origin = [0,Trans.ElementPos(Trans.elemInArr+1,2)*scaleToWvl,0];
    end
    TX(i).Steer = [angRange(i),0];
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(1));
end

lastBmodeTxArr1 = P.numPlanes;



for i = 1:dop.numBeams
    TX(lastBmodeTxArr1+i).waveform = 2;
    [~,centerElem] = min(abs(Trans.ElementPos(1:Trans.elemInArr,1)*scaleToWvl-centerFocus(i)));
    if useArr1
        TX(lastBmodeTxArr1+i).Origin = ([Trans.ElementPos(centerElem,1),Trans.ElementPos(1,2),0])*scaleToWvl;
        elem2Use = (Trans.ElementPos(1:Trans.elemInArr,1)*scaleToWvl > (TX(lastBmodeTxArr1+i).Origin(1)-elemInAp/2*Trans.spacing)) & (Trans.ElementPos(1:Trans.elemInArr,1)*scaleToWvl < TX(lastBmodeTxArr1+i).Origin(1)+elemInAp/2*Trans.spacing);
        TX(lastBmodeTxArr1+i).Apod(1:Trans.elemInArr) = elem2Use;

    else
        TX(lastBmodeTxArr1+i).Origin = [Trans.ElementPos(centerElem,1),Trans.ElementPos(Trans.elemInArr+1,2),0];
        
        elem2Use = ((Trans.ElementPos((Trans.elemInArr+1):(2*Trans.elemInArr),1)*scaleToWvl) > (TX(lastBmodeTxArr1+i).Origin(1)-elemInAp/2*Trans.spacing)) & (Trans.ElementPos((Trans.elemInArr+1):(2*Trans.elemInArr),1)*scaleToWvl < TX(lastBmodeTxArr1+i).Origin(1)+elemInAp/2*Trans.spacing);
        TX(lastBmodeTxArr1+i).Apod((Trans.elemInArr+1):(2*Trans.elemInArr)) = elem2Use;
        
    end
    TX(lastBmodeTxArr1+i).Steer = [0,0];
    TX(lastBmodeTxArr1+i).focus = dop.focalDepth;
    TX(lastBmodeTxArr1+i).Delay = computeTXDelays(TX(lastBmodeTxArr1+i));
    TX(lastBmodeTxArr1+i).peakBLMax = 11;
    TX(lastBmodeTxArr1+i).peakCutOff = 0;
    TX(lastBmodeTxArr1+i).TXPD = computeTXPD(TX(lastBmodeTxArr1+i),PData(2));
end

lastDopTxArr1 = lastBmodeTxArr1 + dop.numBeams;

%% TPC
TPC(1).name = 'US Power';
TPC(1).maxHighVoltage = 25;
TPC(1).hv = 10;

%%
% Specify TGC Waveform structure.
% TGC(1).CntrlPts = [189,314,457,698,770,911,948,976];
TGC(1).CntrlPts = ones(1,8)*1023;
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

TGC(2).CntrlPts = [189,314,457,698,770,911,948,976];
TGC(2).rangeMax = dop.endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));
%% Specify Receive structure arrays.
% - We need m Receive structures for each frame.
maxAcqLength = ceil(sqrt(P.endDepth^2 + P.fovWidth^2));
Receive = repmat(struct('Apod', zeros(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode','BS100BW', ...
                        'mode', 1, ...
                        'callMediaFunc', 0), 1, (Resource.RcvBuffer(1).numFrames*P.numAvgs*P.numPlanes + Resource.RcvBuffer(2).numFrames*dop.numBeams*dop.PRI));
% - Set event specific Receive attributes.
for i = 1:Resource.RcvBuffer(1).numFrames
    for j = 1:P.numPlanes
        for k = 1:P.numAvgs
            if k == 1
                Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).mode = 0;
            end
            if useArr1
%                 Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).Apod(1:Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
                Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).Apod(1:Trans.elemInArr) = ones(1,Trans.elemInArr);
            else
%                 Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).Apod(Trans.elemInArr+1:2*Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
                Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).Apod(Trans.elemInArr+1:2*Trans.elemInArr) = ones(1,Trans.elemInArr);

            end
            Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).framenum = i;
            Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).acqNum = (j-1)*P.numAvgs+k;
        end
    end
end
lastBmodeRcvArr1 = Resource.RcvBuffer(1).numFrames*P.numPlanes*P.numAvgs;

for i = 1:Resource.RcvBuffer(2).numFrames
    for j = 1:dop.numBeams*dop.PRI
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).bufnum = 2;
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).mode = 0;
        if useArr1
            Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).Apod(1:Trans.elemInArr) = ones(1,Trans.elemInArr);
        else
            Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).Apod(Trans.elemInArr+1:2*Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
        end
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).framenum = i;
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).acqNum = j;
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).endDepth = dop.endDepth;
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).startDepth = dop.startDepth;
        load('BP8_4_12_4MHz.mat','BP8_4to12_4MHz');
        Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).InputFilter = BP8_4to12_4MHz(1:21);
%         Receive(lastBmodeRcvArr1+(i-1)*dop.numBeams*dop.PRI+j).sampleMode = 'BS100BW';
    end
end
lastDopRcvArr1 = lastBmodeRcvArr1 + Resource.RcvBuffer(2).numFrames*dop.numBeams*dop.PRI;

% for i = 1:Resource.RcvBuffer(3).numFrames
%     for j = 1:P.numPlanes
%         for k = 1:P.numAvgs
%             if k == 1
%                 Receive((i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).mode = 0;
%             end
%             Receive(lastDopRcvArr1+(i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).Apod(Trans.elemInArr+1:2*Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
%             Receive(lastDopRcvArr1+(i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).framenum = i;
%             Receive(lastDopRcvArr1+(i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).acqNum = (j-1)*P.numAvgs+k;
%             Receive(lastDopRcvArr1+(i-1)*P.numAvgs*P.numPlanes + (j-1)*P.numAvgs + k).bufnum = 3;
%         end
%     end
% end
% lastBmodeRcvArr2 = lastDopRcvArr1+Resource.RcvBuffer(3).numFrames*P.numPlanes*P.numAvgs;
% 
% for i = 1:Resource.RcvBuffer(4).numFrames
%     for j = 1:dop.numAngs*dop.PRI
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).bufnum = 4;
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).mode = 0;
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).Apod(Trans.elemInArr+1:2*Trans.elemInArr) = kaiser(Trans.elemInArr,P.beta);
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).framenum = i;
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).acqNum = j;
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).endDepth = dop.endDepth;
%         Receive(lastBmodeRcvArr2+(i-1)*dop.numAngs*dop.PRI+j).endDepth = dop.startDepth;
%     end
% end
% lastDopRcvArr2 = lastBmodeRcvArr2 + Resource.RcvBuffer(2).numFrames*dop.numAngs*dop.PRI;
%% RcvProfile
RcvProfile(1).AntiAliasCutoff = 20;
RcvProfile(1).PgaGain = 30;
RcvProfile(1).LnaGain = 24;
RcvProfile(1).LnaHPF = 200;
RcvProfile(1).LnaZinSel = 2;

RcvProfile(2).AntiAliasCutoff = 20;
RcvProfile(2).PgaGain = 24; %30
RcvProfile(2).LnaGain = 24;
RcvProfile(2).LnaHPF = 200;
RcvProfile(2).LnaZinSel = 30;

%% Specify Recon structure arrays.
% - We need one Recon structure.  Each frame will use
%   m ReconInfo structures, since we are using m
%   synthetic aperture acquisitions.
Recon = repmat(struct('senscutoff', 0.6, ...
               'pdatanum', 1, ...
               'rcvBufFrame', -1, ...     % use most recently transferred frame
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...  % auto-increment ImageBuffer each recon
               'RINums', []),1,2);

%Array 1 Doppler
Recon(2).pdatanum = 2;
Recon(2).IntBufDest = [2,-1];
Recon(2).ImgBufDest = [0, 0];
% Recon(2).newFrameTimeout = 10000;

% %Array 2 Bmode
% Recon(3).IntBufDest = [3,1];
% Recon(3).ImgBufDest = [3,-1];
% 
% %Array 2 Doppler
% Recon(4).IntBufDest = [4,-1];
% Recon(4).ImgBufDest = [4,-1];


%% Define ReconInfo structures.
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % accumulate IQ data.
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'regionnum', 1), 1, (P.numPlanes+dop.numBeams*dop.PRI));
% - Set specific ReconInfo attributes.

ReconInfo(1).Pre = 'clearInterBuf';
ReconInfo(1).mode = 'replaceIQ';  % replace IQ data
for j = 1:P.numPlanes  % For each row in the column
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
% ReconInfo(P.numPlanes).mode = 'accumIQ_replaceIntensity'; % accum & detect
ReconInfo(P.numPlanes).Post = 'IQ2IntensityImageBuf'; % accum & detect
lastBmodeRIArr1 = P.numPlanes;
Recon(1).RINums = 1:lastBmodeRIArr1;

%Doppler Arr 1 RI
ReconInfo(lastBmodeRIArr1+1).Pre = 'clearInterBuf';
ReconInfo(lastBmodeRIArr1+1).mode = 'replaceIQ';
for k = 1:dop.PRI
    for j = 1:dop.numBeams  % For each row in the column
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numBeams+j).txnum = lastBmodeTxArr1+j;
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numBeams+j).rcvnum = lastBmodeRcvArr1+(k-1)*dop.numBeams+j;
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numBeams+j).pagenum = k;
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numBeams+j).mode = 'replaceIQ';
        ReconInfo(lastBmodeRIArr1+(k-1)*dop.numBeams+j).regionnum = j;
    end
end
lastDopRIArr1 = lastBmodeRIArr1 + dop.numBeams*dop.PRI;
Recon(2).RINums = lastBmodeRIArr1+1:lastDopRIArr1;

% % Bmode Arr 2
% ReconInfo(lastDopRIArr1+1).Pre = 'clearInterBuf';
% ReconInfo(lastDopRIArr1+1).mode = 'replaceIQ';  % replace IQ data
% for j = 1:P.numPlanes  % For each row in the column
%     ReconInfo(lastDopRIArr1+j).txnum = lastDopTxArr1+j;
%     ReconInfo(lastDopRIArr1+j).rcvnum = lastDopRcvArr1+j;
%     ReconInfo(lastDopRIArr1+j).regionnum = j;
% end
% % ReconInfo(P.numPlanes).mode = 'accumIQ_replaceIntensity'; % accum & detect
% ReconInfo(lastDopRIArr1+P.numPlanes).Post = 'IQ2IntensityImageBuf'; % accum & detect
% lastBmodeRIArr2 = lastDopRIArr1+P.numPlanes;
% Recon(3).RINums = lastDopRIArr1+1:lastBmodeRIArr2;
% 
% %Doppler Arr 2 RI
% ReconInfo(lastBmodeRIArr2+1).Pre = 'clearInterBuf';
% ReconInfo(lastBmodeRIArr2+1).mode = 'replaceIQ';
% for k = 1:dop.PRI
%     for j = 1:dop.numAngs  % For each row in the column
%         ReconInfo(lastBmodeRIArr2+(k-1)*dop.numAngs+j).txnum = lastBmodeTxArr2+j;
%         ReconInfo(lastBmodeRIArr2+(k-1)*dop.numAngs+j).rcvnum = lastBmodeRcvArr2+(k-1)*dop.numAngs+j;
%         ReconInfo(lastBmodeRIArr2+(k-1)*dop.numAngs+j).regionnum = j;
%         ReconInfo(lastBmodeRIArr2+(k-1)*dop.numAngs+j).pagenum = k;
%         ReconInfo(lastBmodeRIArr2+(k-1)*dop.numAngs+j).mode = 'replaceIQ';
%     end
% end
% lastDopRIArr2 = lastBmodeRIArr2 + dop.numAngs*dop.PRI;
% Recon(4).RINums = lastBmodeRIArr2+1:lastDopRIArr2;

%% Specify Process structure array.
pers = 20;
dopThresh = 0.01;
maxPwr = 20;
compFactorB = 40;
rejB = 0;
pGainB = 2;

compFactorD = 70;
rejD = 2;
pGainD = 50;
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

CDIProcess = 2;
Process(2).classname = 'Doppler';
Process(2).method = 'computeCFIFreqEst';
Process(2).Parameters = {'IntBufSrc',[2,-1],...   
                         'SrcPages',[1 dop.PRI-1],...   
                         'IntBufDest',[3,1],...    
                         'ImgBufDest',[2,-1],...            
                         'pdatanum',2,...      
                         'prf',ceil(dop.PRF),...
                         'wallFilter','none',...
                         'pwrThreshold',dopThresh,...
                         'maxPower',maxPwr};

bkgdBmodeProcess = 3;
dispThresh = 5;
%Display Bkgd Ultrasound Arr 1
Process(3).classname = 'Image';
Process(3).method = 'imageDisplay';
Process(3).Parameters = {'imgbufnum',1,...   % number of buffer to process.
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
                         'displayWindow',2,...
                         'threshold',dispThresh};

%Display Doppler Arr 1
Process(4).classname = 'Image';
Process(4).method = 'imageDisplay';
Process(4).Parameters = {'imgbufnum',2,...   % number of buffer to process.
                         'framenum', -1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',pGainD,...            % pgain is image processing gain
                         'reject',rejD,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',compFactorD,...
                         'mappingMethod','upperHalf',...
                         'display',1,...      % display image after processing
                         'displayWindow',2,...
                         'srcData','signedColor'};

Process(5).classname = 'External';
Process(5).method = 'myTic';
Process(5).Parameters = {'srcbuffer','none',...
                         'dstbuffer','none'};
                     
                     
PDIProcess = 6;
Process(6).classname = 'Doppler';
Process(6).method = 'computeCFIPowerEst';
Process(6).Parameters = {'IntBufSrc',[2,-1],...   
                         'SrcPages',[1 dop.PRI-1],...   
                         'ImgBufDest',[3,-1],...            
                         'pdatanum',2,...      
                         'prf',ceil(dop.PRF),...
                         'wallFilter','none',...
                         'pwrThreshold',dopThresh,...
                         'maxPower',maxPwr};
                     
%Display Pwr Doppler Arr 1
Process(7).classname = 'Image';
Process(7).method = 'imageDisplay';
Process(7).Parameters = {'imgbufnum',3,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',2,...    % number of PData structure to use
                         'pgain',pGainD,...            % pgain is image processing gain
                         'reject',rejD,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',compFactorD,...
                         'mappingMethod','upperHalf',...
                         'display',1,...      % display image after processing
                         'displayWindow',3,...
                         'srcData','unsignedColor'};
                     
%Display Bkgd Ultrasound Arr 1
bkgdPwrDopBmodeProcess = 8;
Process(8).classname = 'Image';
Process(8).method = 'imageDisplay';
Process(8).Parameters = {'imgbufnum',1,...   % number of buffer to process.
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
                         'displayWindow',3,...
                         'threshold',dispThresh};
                     
% %Display Bkgd Ultrasound Arr 2
% Process(5).classname = 'Image';
% Process(5).method = 'imageDisplay';
% Process(5).Parameters = {'imgbufnum',3,...   % number of buffer to process.
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
% Process(6).classname = 'Image';
% Process(6).method = 'imageDisplay';
% Process(6).Parameters = {'imgbufnum',4,...   % number of buffer to process.
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
%% Specify SeqControl structure arrays.
SeqControl(1).command = 'jump'; % jump back to start.
SeqControl(1).argument = 1;

SeqControl(2).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(2).argument = 1/(P.maxFPS*P.numPlanes)*10^6;  % 200 usec
SeqControl(3).command = 'returnToMatlab';

SeqControl(4).command = 'sync';
SeqControl(4).argument = 5000000;

SeqControl(5).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(5).argument = 50000;  % 200 usec

sc_dop_ttna = 6;
SeqControl(6).command = 'timeToNextAcq';  % time between synthetic aperture acquisitions
SeqControl(6).argument = ceil(1/(dop.PRF*dop.numBeams)*1E6);

SeqControl(7).command = 'setRcvProfile';  % time between synthetic aperture acquisitions
SeqControl(7).argument = 1; 

SeqControl(8).command = 'setTPCProfile';  % time between synthetic aperture acquisitions
SeqControl(8).argument = 1;
SeqControl(8).condition = 'immediate';

SeqControl(9).command = 'noop';  % time between synthetic aperture acquisitions
SeqControl(9).argument = 100000*0.2;

SeqControl(10).command = 'setRcvProfile';  % time between synthetic aperture acquisitions
SeqControl(10).argument = 2; 

SeqControl(11).command = 'jump'; % jump back to start.

nsc = 12; % nsc is count of SeqControl objects

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

    Event(n).info = 'Store PDI Bkgd US';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 3;
    Event(n).seqControl = 0;
    n = n+1;
    
    Event(n).info = 'Store PDI Bkgd US ';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 0;
    Event(n).process = 8;
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
Event(n).info = 'Set Doppler Rcv Profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 10;
n = n+1;

for i = 1:dop.PRI
    for j = 1:dop.numBeams
        Event(n).info = 'Doppler Plane Wave TX';
        Event(n).tx = lastBmodeTxArr1 + j;
        Event(n).rcv = lastBmodeRcvArr1 + (i-1)*dop.numBeams + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 6;
        n = n+1;
    end
end
Event(n-1).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    
Event(n).info = 'Recon Doppler';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 2;
Event(n).process = 2;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Display Doppler';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 4;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

startPwrDop = n;
Event(n).info = 'Set Doppler Rcv Profile';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 10;
n = n+1;

for i = 1:dop.PRI
    for j = 1:dop.numBeams
        Event(n).info = 'Doppler Focused Wave TX';
        Event(n).tx = lastBmodeTxArr1 + j;
        Event(n).rcv = lastBmodeRcvArr1 + (i-1)*dop.numBeams + j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = 6;
        n = n+1;
    end
end
Event(n-1).seqControl = nsc;
    SeqControl(nsc).command = 'transferToHost';
    nsc = nsc + 1;
    
Event(n).info = 'Estimate Power Sig';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 2;
Event(n).process = 6;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Display Pwr Doppler';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 7;
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Jump';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = 1;
n = n+1;

%% UI

UI(1).Control = {'UserB6','Style','VsPushButton','Label','CDI'};
UI(1).Callback = text2cell('%SingleAcqDoppler');

UI(2).Control =  {'UserB5','Style','VsSlider','Label','CDI Thresh',...
                  'SliderMinMaxVal',[0,1,dopThresh],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%1.2f'};
UI(2).Callback = text2cell('%changePwrThresh');

UI(3).Control =  {'UserB4','Style','VsSlider','Label','CDI Max Pwr',...
                  'SliderMinMaxVal',[0,100,maxPwr],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%d'};
UI(3).Callback = text2cell('%changeMaxPwr');

UI(4).Control =  {'UserB3','Style','VsSlider','Label','PRF (kHz)',...
                  'SliderMinMaxVal',[0,10,dop.PRF/1000],...
                  'SliderStep',[0.1,1],'ValueFormat','%0.1f'};
UI(4).Callback = text2cell('%changePRF');

UI(5).Control =  {'UserB2','Style','VsSlider','Label','BMode Reject',...
                  'SliderMinMaxVal',[0,100,rejB],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%2d'};
UI(5).Callback = text2cell('%bmodereject');

UI(6).Control =  {'UserC1','Style','VsSlider','Label','BMode Compress',...
                  'SliderMinMaxVal',[0,1,compFactorB/100],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%0.2f'};
UI(6).Callback = text2cell('%changeCompression');

UI(7).Control =  {'UserB1','Style','VsSlider','Label','BMode Gain',...
                  'SliderMinMaxVal',[0,100,pGainB],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%2d'};
UI(7).Callback = text2cell('%changeGain');

UI(8).Control =  {'UserC8','Style','VsSlider','Label','Doppler Focal Depth',...
                  'SliderMinMaxVal',[1,20,dop.focalDepthMm],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%0.2f'};
UI(8).Callback = text2cell('%changeFocalDepth');
UI(9).Control =  {'UserC7','Style','VsSlider','Label','Disp Thresh',...
                  'SliderMinMaxVal',[1,101,dispThresh],...
                  'SliderStep',[0.01,0.1],'ValueFormat','%2d'};
UI(9).Callback = text2cell('%changeDispThresh');

UI(10).Control = {'UserC2','Style','VsPushButton','Label','PDI'};
UI(10).Callback = text2cell('%PwrDopAcq');

% Save all the structures to a .mat file.
filename = 'C:\Users\Administrator\Documents\Vantage-4.5.3-2107301223\MatFiles\MUSIC_FocusedBeamDoppler.mat';
save(filename);
VSX;

% txt1 = input('Filename for saved data?','s');
% txt2 = input('What board was connected to connector L?','s');
% savefast(['20230322_PigSurgery\',txt1,'_Board',txt2,'.mat'],'Event','IData','QData','Receive','Trans','TX','Resource','PData','RcvProfile','Recon','ReconInfo','ImgData','ImgDataP','RcvData')
return

% **** Callback routines to be converted by text2cell function. ****
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

%RangeChangeCallback - Range change
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.endDepth'));
    return
end
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
assignin('base','P',P);
evalin('base','PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));');
evalin('base','PData(1).Region = computeRegions(PData(1));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback
%% UI(2) - Single Acquisition Doppler
%SingleAcqDoppler
singleshot_start = evalin('base','startDopArr1');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',singleshot_start};
assignin('base','Control', Control);

%SingleAcqDoppler

%% UI(3) - change power threshold
%changePwrThresh
pwrThresh = evalin('base','dopThresh');
CDIProcess = evalin('base','CDIProcess');
PDIProcess = evalin('base','PDIProcess');


Process = evalin('base','Process');
for k = 1:2:length(Process(CDIProcess).Parameters)
    if strcmp(Process(CDIProcess).Parameters{k},'pwrThreshold'), Process(CDIProcess).Parameters{k+1} = UIValue; end
end
for k = 1:2:length(Process(PDIProcess).Parameters)
    if strcmp(Process(PDIProcess).Parameters{k},'pwrThreshold'), Process(PDIProcess).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Doppler threshold.
Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',CDIProcess,'pwrThreshold',UIValue};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',PDIProcess,'pwrThreshold',UIValue};

% Control = evalin('base','Control');
% Control(1).Command = 'set&Run';
% Control(1).Parameters = {'Process',dopProcess,'pwrThreshold',UIValue};

pwrThresh = UIValue;
assignin('base','pwrThresh',pwrThresh);
assignin('base','Control', Control);

%changePwrThresh

%% UI(4) - change max power
%changeMaxPwr
maxPwr = evalin('base','maxPwr');
CDIProcess = evalin('base','CDIProcess');
PDIProcess = evalin('base','PDIProcess');

Process = evalin('base','Process');
for k = 1:2:length(Process(CDIProcess).Parameters)
    if strcmp(Process(CDIProcess).Parameters{k},'maxPwr'), Process(CDIProcess).Parameters{k+1} = UIValue; end

end
for k = 1:2:length(Process(PDIProcess).Parameters)
    if strcmp(Process(PDIProcess).Parameters{k},'maxPwr'), Process(PDIProcess).Parameters{k+1} = UIValue; end
end


Control = evalin('base','Control');
Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',CDIProcess,'maxPower',UIValue};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Process',PDIProcess,'maxPower',UIValue};
maxPwr = UIValue;
assignin('base','maxPwr',maxPwr);
assignin('base','Control', Control);

%changeMaxPwr

%% UI(5) - Change PRF
%changePRF
dop = evalin('base','dop');
CDIProcess = evalin('base','CDIProcess');
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
Control(1).Parameters = {'Process',CDIProcess,'prf',dop.PRF};

setLoopCntTo1 = evalin('base','setLoopCntTo1');
SeqControl(sc_dop_ttna).argument = 1/(dop.numBeams*dop.PRF*1E-6);
assignin('base','SeqControl',SeqControl);
Control(2).Command = 'update&Run';
Control(2).Parameters = {'SeqControl'};
Control(3).Command = 'set&Run';
Control(3).Parameters = {'Parameters',1,'startEvent',setLoopCntTo1};

disp([num2str(dop.PRI*((dop.numBeams)*SeqControl(sc_dop_ttna).argument)/1000),...
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

%changeFocalDepth
dop = evalin('base','dop');
CDIProcess = evalin('base','CDIProcess');
TX = evalin('base','TX');
Trans = evalin('base','Trans');
centerFocus = evalin('base','centerFocus');
PData = evalin('base','PData');
scaleToWvl = 1/Trans.wvl;
lastBmodeTxArr1 = evalin('base','lastBmodeTxArr1');
dop.focalDepthMm = UIValue;
dop.focalDepth = dop.focalDepthMm/Trans.wvl;
startBMode = evalin('base','startBMode');

elemInAp = ceil((dop.focalDepth/dop.fNum)/Trans.spacing);

if elemInAp > Trans.elemInArr
    elemInAp = Trans.elemInArr;
end


for i = lastBmodeTxArr1+1:lastBmodeTxArr1+dop.numBeams
    elem = i-lastBmodeTxArr1;
    [~,centerElem] = min(abs(Trans.ElementPos(1:Trans.elemInArr,1)*scaleToWvl-centerFocus(elem)));

    TX(i).Origin = ([Trans.ElementPos(centerElem,1),Trans.ElementPos(1,2),0])*scaleToWvl;
    elem2Use = (Trans.ElementPos(1:Trans.elemInArr,1)*scaleToWvl > (TX(i).Origin(1)-elemInAp/2*Trans.spacing)) & (Trans.ElementPos(1:Trans.elemInArr,1)*scaleToWvl < TX(i).Origin(1)+elemInAp/2*Trans.spacing);
    TX(i).Apod(1:Trans.elemInArr) = elem2Use;
    TX(i).focus = dop.focalDepthMm/Trans.wvl;
    TX(i).Delay = computeTXDelays(TX(i));
    TX(i).TXPD = computeTXPD(TX(i),PData(2));
end

assignin('base','TX',TX);
assignin('base','dop',dop);

Control = evalin('base','Control');
Control(1).Command = 'update&Run';
Control(1).Parameters = {'TX'};
Control(2).Command = 'set&Run';
Control(2).Parameters = {'Parameters',1,'startEvent',startBMode};
assignin('base','Control',Control);

%changeFocalDepth

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

%PwrDopAcq
startPwrDop = evalin('base','startPwrDop');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',startPwrDop};
assignin('base','Control', Control);
%PwrDopAcq