clear all

%% BMode User Parameters
P.startDepthMm = 0;   % start depth in mm
P.endDepthMm = 25;    % end depth in mm
P.fovWidthMm = 15;    % the width of the bottom of the image (mm)
P.numAvgs = 1;        % number of averages for the Bmode image
P.numPlanes = 10;    % number of planes used to create a Bmode image
P.maxFPS = 100;        % Frame per second for Bmode imaging
P.c = 1540;           % Speed of Sound (m/s)      
P.beta = 4;         % Narrowness of the kaiser window
P.useElem = [1, 64];
%%
fus.totaldur = 120; %seconds %No use
fus.PRF = 1000; %Hz
fus.burst = .35 %ms
fus.V = 2.5; %volts
fus.freq = 2.5E6;
fus.focalMm = [0,0,7]; %mm
fus.useElem = [209 224];%[209 224], [193 208] , [193 224] Change lines 180 and 185

%% Specify system parameters.
Resource.Parameters.numTransmit = 244;     % number of transmit channels.
Resource.Parameters.numRcvChannels = 244;  % number of receive channels.
Resource.Parameters.speedOfSound = P.c;
Resource.Parameters.verbose = 2;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.waitForProcessing = 0; % DMA transfers wait for processing.
Resource.Parameters.simulateMode = 0;
Resource.Parameters.Connector = [1,2];

% HIFU Codes
Resource.HIFU.externalHifuPwr = 1; % Use External Power
Resource.HIFU.extPwrComPortID = 'COM5'; % Subject to changes
Resource.HIFU.psType = 'QPX600DP'; % set to 'QPX600DP' to match supply being used
%
%% Specify Trans structure array.
Trans.name = 'MUSIC3_Biobox';
Trans.units = 'wavelengths'; % required in Gen3 to prevent default to mm units
Trans.frequency = 2.5;
TWFreq = 250/ceil(250/(4*Trans.frequency))/4;
Trans.frequency = TWFreq;
Trans.wvl = Resource.Parameters.speedOfSound*1000/(Trans.frequency*10^6);
Trans.Bandwidth = [7 14];
Trans = computeTrans2(Trans);  % L11-4v transducer is 'known' transducer so we can use computeTrans.
Trans.maxHighVoltage = 95;  % set maximum high voltage limit for pulser supply.
Trans.name = 'custom';

%
Trans.impedance = 100; %Set for HIFU config
%

%%
fus.focal = fus.focalMm/((P.c/fus.freq)*1000)+[0,mean(Trans.ElementPos(fus.useElem(1):fus.useElem(2),2)),0];
fus.totalPulses = fus.totaldur*fus.PRF;
fus.maxLoopCnt = 65535;

%%
% TW(1).type = 'envelope';
% assert(ceil(fus.burst*(fus.freq)/1000) < 10000, 'Burst length is too long')

% TW(1).envNumCycles = ceil(fus.burst*(fus.freq)/1000)/2; %10;
% TW(1).envFrequency = ones(1,TW(1).envNumCycles)*fus.freq/1E6;
% TW(1).envPulseWidth = ones(1,TW(1).envNumCycles);

TW(1).type = 'parametric';
TW(1).Parameters = [2.4, 1,2,1];