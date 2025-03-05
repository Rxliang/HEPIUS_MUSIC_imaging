% Code to generate ultrasound sequence, and implement some UI
% callback functions

systemID = getenv('BNS_SYSTEM_ID');

if ~strcmp(systemID, 'veraBCM')
 % deal with exception at Baylor
 %setsystemid256
 setsystemid64
end

clear
close all
clear all
close all
clear all
clear saveiq sdlargefn sdlargeoutfn

% Default version info, overwritten by productionsettings.m
version = [];
version.name = 'BNS205 *** Not for clinical use ***';
version.number = '0.0';

simulateMode = 0; % hardware provides realtime data
%simulateMode = 1; % Verasonics default data (no relevance to this application)
%simulateMode = 2; % saved real data (can be slow to load)

debugSeqSwitchMode=0;

mfile = mfilename;

% check we don't use unsupported code
execState = 'production';
execState = 'walltrack';
%execState = 'field measurement';
%execState = 'phantom';

productionFlag = strcmp(execState, 'production');
fieldMeasurementFlag = strcmp(execState, 'field measurement');
phantomFlag = strcmp(execState, 'phantom');

% automatically run VSX to launch sequence after script has completed
autoVSX=0;
% don't hide VSX: usgui will override this and hide it when it is run
Mcr_GuiHide = 0;

% setup paths

% If environment variables are not defined for system,
% for development only, you can run this to define them
if 0
  uscmsg(execState);
  definebnsenv
end

veraPath=getenv('VERASONICS_VPF_ROOT');
if isempty(veraPath)
  activate
  veraPath=getenv('VERASONICS_VPF_ROOT');
end

vtHome = getenv('MUSIC_HOME');
vtData = getenv('MUSIC_DATA');

checkbnspaths;

if ~exist(vtHome, 'dir')
  mkdir(vtHome);
end

if ~exist(vtData, 'dir')
  mkdir(vtData);
end

tmpDataPath = [vtData '/tmp'];
if ~exist(tmpDataPath, 'dir')
  mkdir(tmpDataPath);
end

tmpDirSub = fullfile('nicp','tmp');
tmpDir = fullfile(vtData, tmpDirSub);
if ~exist(tmpDir, 'dir')
  mkdir(tmpDir);
end

% storage path for IQ data
iqDataPath = [vtData '/nicp/iqdata/'];
if ~exist(iqDataPath)
  mkdir(iqDataPath)
end

matOutPath = [veraPath];

dateStr = datestr(now, 'yyyymmdd_HHMM');

% Some UI parameters
evParms = [];
evParms.lock.ctrlLockSleepTime_s = 0.25;
evParms.lock.ctrlUpdate = 0;
evParms.flag.debugSeqSwitchMode=debugSeqSwitchMode;
evParms.UI.waitbarHandle = 0;

% set known elements of Resource structure

Resource.Parameters.speedOfSound = 1540;
Resource.Parameters.speedCorrectionFactor = 1.0;
Resource.Parameters.initializeOnly = 0;
Resource.Parameters.simulateMode = simulateMode;

matFileSuffixSim = '';

if Resource.Parameters.simulateMode 
  matFileSuffixSim = ['_sim' num2str(Resource.Parameters.simulateMode)];
end
               
%  Resource.Parameters.simulateMode = 1 %forces simulate mode, even if hardware is present.
%  Resource.Parameters.simulateMode = 2 %stops sequence and processes RcvData continuously.

if 0
  uscmsg(execState);
  % This is to support simulation mode
  % Specify Media object.
  pt1;
  Media.function = 'movePoints';
end

% Debugging tool section

% Vantage software has memory leak that causes crash for large
% number of CDI PRIs. Verasonics has not fixed this, so we work 
% within the limit. This code allows demonstration of this bug
% and may be used to verify any eventual resolution of the issue

demoWallFilterIssue=0;
if demoWallFilterIssue
  disp('Warning: CDI PRI issue debug setting');
  % Cases to illustrate issues with wall filter with large number of PRIs
  demoStableCase1 = 1; % 77 src pages
  demoStableCase2 = 0; % 98 src pages, no wall filter
  demoUnstableCase = 0; % 78 src pages, with wall filter
  % allow code to run on a generic PC
  setenv('BNS_HOME', 'c:\temp');
  setenv('BNS_DATA', 'c:\temp');
else
  demoStableCase1 = 0;
  demoStableCase2 = 0;
  demoUnstableCase = 0;
end

% allow debugging of aspects of this script by inserting pauses
useCheckpoints=0;
% on systems with a Trigger license, enable debugging using
% hardware trigger signals. Also needed for ultrasound field measurements
sendTriggers=0; % also sync at start of acq
% trigger specification
trig=[];
trig.frameNos = 0; % frame numbers to trigger on 
trig.onFirstPulseOfTypeOnly = 0;
trig.onB = 0;
trig.onSD= 0;
trig.onCFI= 0;
trig.onBPreCFI=0; % on the B mode pulses before CFI seq

% For debugging, can force recording of spectral Doppler IQ
forceRecord=0; % don't require semaphore to record

% B mode is generated using different transmission pulses in
% non-SDOnlyInPlace mode
separateBAcq=1;

% transducer specification
% default is to use GE M5Sc-D with 64 channel MUX system
usem5 = 1; % Use GE M5Sc-D
usem5_64 = 1; % use 64 central channels of central row of M5Sc-D
              % using ATL connector and adapter
% these are not supported and are not tested
use6s = 0; % use GE 6S-D
use312 = 0; % use GE L3-12D
usec15 = 0; % use GE C1-5D
usec16 = 0; % use GE C1-5D
if usec16
  usec15=1;
end

% default system uses UTA260MUX front panel and adapter for GE probe
UTA260MUX=2; % set to 2 for Vantage 4.0.0 and later

Trans = [];
% default frequency, used for Doppler

Trans.frequency = 2.8409; 
%Trans.frequency = 5.6818; 

% B-mode parameters

%BModeFrequency = 4.6;
BModeFrequency = Trans.frequency;
%BModeFrequency = 6.9444;
%BModeFrequency = 5.6818; % also uses 22 AD
%BModeFrequency = Trans.frequency; % also uses 22 ADC

% long ttna at end of last non cfi event in frame
evParms.ev.lastFrameTTNA = 0;

evParms.ev.transferTTNA_us = 10000;

% after last Doppler tx when regular CFI not switched on

evParms.ev.transferTTNA_postSD_us = 10000; % 60000;

%BModeFrequency = 3.8;
BSteerAngle_deg = 0; % worth making adjustable in future  

% allows spectral Doppler signals to be added together along ray
% (simulating a CW-like measurement). Not supported and not tested.
narrowbanding=0;
SDOnly=0; % do not create sequence for SDOnly mode (no B, no CFI)
duplexMode=1; % B+SD simultaneously, which can be switched to B+CDI
numSD=2;

evParms.B.depthOffset_mm = 0;
evParms.state.initDepthSliderDone = 0;

% spectrogram parameters. these are read from base workspace by
% spectralDoppler scripts using evalin
evParms.SD.functionName = 'sdwrapfn';
evParms.SD.enableInstances = [];
evParms.SD.forceInit = [0 0]; % ch 1,2 used to force reinitialization
evParms.SD.despeckle = 0.4;
evParms.SD.dopSweepTimeIndex = 2; %1-4 avail for 20msec frame time
evParms.SD.noisePersist=0.95;
evParms.SD.SDDynRangeDB = 13;% not including FFT processing gain
evParms.SD.SDDynRangeDB = 1;% not including FFT processing gain
evParms.SD.SDDynRangeDB = 0.25;% not including FFT processing gain
%baseLineShiftNorm = -0.375;
evParms.SD.baseLineShiftNorm = -0.2;
% default wall filter cutoff, normalized
 % wall filter normalized cutoff frequency
evParms.default.SD.dopWFCutoffNorm = 0.025*2;
evParms.SD.dopWFCutoffNorm = evParms.default.SD.dopWFCutoffNorm;
evParms.SD.titleSuffix = ' marker';
evParms.SD.markerColorPower = [51/255   51/255  1; 
                               0.2   1 .2];
evParms.SD.markerColorFlow = [ 1   1   0  ; 
                               0.2   1 .2];

% use a focussed beam (1) or a pure plane wave (0).
focussedDoppler=1;
if focussedDoppler 
  % optimized to get narrow beam at > 6 cm using GE M5Sc-D with 64 channels
  txFocus_mm = 160; 
end

% IQ data saving parameters
numSaveFramesPerBatch = 2000; %650; % 100 is 5 sec
SDLargeFramesPerFile = 500;
forceHighPRFModeForRecordingIQ = 1; % productionsettings will override

% ultrasound parameters

% sequence timing. 
numFrames=10; % per sequence
if Resource.Parameters.simulateMode==2 
  numFrames = 100; % load/save about 5 seconds of data
end

framePeriod_s = 50e-3; % 20 frames/sec
evParms.ev.watchdogInterval_s = 2;
% make sure that the PCC is updated with consisent with this
% timeout value
evParms.ev.watchdogInterval_frames = ...
    round(evParms.ev.watchdogInterval_s/framePeriod_s);
evParms.ev.watchdogFilename = fullfile(tmpDirSub, 'watchdogTimestampFile.txt');

% We want a spectral Doppler PRF of 4500 Hz to capture maximum
% velocities encountered in the OA. Since this affects size of
% PData for SD, make as small as possible to improve performance
maxPRFNominal_Hz = 4525; % max PRF for B+SD 
%maxPRFNominal_Hz = 2500; % max PRF for B+SD 
% must be even to accommodate B/SD
maxPRIs = 2*floor(framePeriod_s*maxPRFNominal_Hz/2);
% actual maximum PRF
maxPRF_Hz = maxPRIs/framePeriod_s;
% transit time and event overhead limit this. For a maximum depth
% of 90cm we can tolerate only a limited duplex
% PRF. Experimentally, this is an acceptable delay between B and SD
% pulses for this depth. Corresponds to about 3.3 kHz.
PRIAcq_uSec = 150.5; % sequence TTNA
% however, at present, we wish to keep number of events constant,
% so we choose a lower PRF for duplex that gives 4.5kHz in
% SDOnlyInPlace mode. This mode provides low quality B-mode imaging
% using SD PRIs, but double the PRF
%PRIAcq_uSec = 221; % 4.5kHz with one wave type
%PRIAcq_uSec = round(1e6/maxPRF_Hz)
%PRIAcqInPlace_uSec = 221;
% transmit waveform parameters
halfCyclesSD=14; % half cycles for SD
%halfCyclesSD=4; % half cycles for SD
halfCyclesB=2; % half cycles for B
halfCyclesCFI=8; 

% enable color Doppler overlay that uses the existing SD PRIs
largeSD=1; 

% start depth for all ultrasound reconstructions
startDepth_mm = 0;

% default spectral Doppler settings
% lateral width of sample volumes
xSVWidth_mm = 1;
% axial length of sample volumes
zSVWidth_mm = 3;
% starting positions of sample volumes [
xSVStart_mm = [-1.237 2.2751];
zSVStart_mm = [56.5381 42.9481];
SDGateDisplayShape = 'rectangle';
% for saving IQ:
xSVWidthIQ_mm = 1;
zSVWidthIQ_mm = 15;

% position of beam center of SD beam
xDopplerBeam_mm=0; % zero is central element

% SD overlay parameters
SDLarge = [];
% default start depth of overlay
SDLarge.zSDStart_mm = 40; % note, used for pdata even when no
                          % sdlarge processing (largeSD=0)

% default lateral position of overlay center
SDLarge.xSDOffset_mm = 0;
% default start depth of overlay
SDLarge.zSDStart_mm = 35;

% These apply in Imaging priority mode:
% default lateral width of overlay
SDLarge.xSDWidth_mm = 25;
% default height of overlay
SDLarge.zSDHeight_mm = 45; 

SDLarge.zSDEnd_mm = SDLarge.zSDStart_mm+SDLarge.zSDHeight_mm;

% For spectrogram mode, we need to ensure we don't have glitching,
% so we limit the area of the overlay

SDLarge.maxSpgArea_mm2 = 25*45;

% save overlay IQ data to file
largeSDParms.largeSDSave=0;
% gen and display largeSD every N frame, set > numFrames to disable
largeSDParms.largeSDFrameInterval=2; %21;
largeSDParms.largeSDFrames = 1:largeSDParms.largeSDFrameInterval:numFrames;
% default CDI mode: power Doppler. *** Flow not implemented yet. ***
largeSDParms.SDLargeCDIMode =  'power';
% default overlay user CDI parameters
largeSDParms.startPRI = 3; % first SD PRI to use
largeSDParms.endPRI = 79; % last SD PRI to use, cannot be > 79 due
                          % to Verasonics bug
largeSDParms.pwrThres = 0.035; % default power threshold
largeSDParms.cpl = 60; % default color priority level
largeSDParms.maxPower=10; % default power to override B-mode
largeSDParms.perspwr=40; % default power Doppler persistence
largeSDParms.SDLargeDisplayOverlay=1; % include sequence event to display the overlay
largeSDParms.wallFilter = 'FIRHigh'; % wall filter that is used for
                                     % overlay only

% SD overlay limits for UI
largeSDParms.maxPowerThreshLimits = [10 100];
largeSDParms.pwrThreshLimits = [0 0.1];
largeSDParms.perspwrLimits = [0 99];
largeSDParms.cplLimits = [0 255];  
SDLarge.zSDMinHeight_mm = 10;  
SDLarge.xSDMinWidth_mm = 10;    
  
% overlay debug parameters
largeSDParms.doSDLargeInternalProc=1; % use Verasonics optimized code to
                                      % calculate overlay, 0
                                      % setting uses Matlab code
                                      % and is not for production use
largeSDParms.procSDLargeOut=0; % save SD overlay image to
                               % intermediate image buffer, not for production
% overlay internal parameters
largeSDParms.priSkip=2; % always use the duplex PRF, or half the
                        % SDOnlyInPlace PRF for the overlay

% CFI mode parameters. These are not used for duplex mode or SD
% overlay
numCFIPRIs = 16; % number of pulses per frame for CDI
dopPRFCFI = 2.5e+03;
pwrThresCFI = 0.35; % jsm, 0.42; % orig
persfreqCFI = 20;
perspwrCFI = 20;
cplCFI = 60;  % color priority level at which 2D overwrites Doppler

regionOverlapFactor3ray=1.1; % PData region overlap when 3 wide
                             % beams are used
regionOverlapFactor4ray=1.2; % PData region overlap when 4 wide
                             % beams are used
regionExtensionFactor=1.4; % Overall PData mag factor vs
                           % txd width


% legacy mode switches, fixed at default values 
startInBCFI = 0; % set true to start with CFI running
doCFI=1; % Enables construction of CFI sequence elements
         % These do not need to be used.
doCFISep = 1; % make infrastructure for CFI to run as separate
              % sequence. it is very hard to run it within a
              % triplex mode, so this is only supported mode
doCFISepAtStart = 0; % set true to start in B/CFI mode with no SD
if doCFISepAtStart
  duplexMode=0;
  SDOnly=0;
end

% development mode switches
multiAngleB=0;
chirpOn=0; % for SD only now
shallowDepth=0; % useful for shirking RF size
% delay HW before each frame
preFrameHWDelay_s = 0;
syncFrameStart = preFrameHWDelay_s; % also sync at start of acq
syncFrameStart = 0;
offlineSDProcFlag = 0; % use spectralDoppler script for offline
                   % processing. Should always be zero.


if 0 % settings for clinical 201807
  Trans.maxHighVoltage = 15;    % set a reasonable high voltage limit.
  UTA260MUX=0;  
  PRIAcq_uSec = 146; % sequence TTNA
  numFrames=10; % per sequence
  numSaveFramesPerBatch=100;
  Trans.frequency = 2.8409;
  startDepth_mm = 0.2;
  xDopplerBeam_mm=0; % zero is central element
  xSVWidth_mm = 1;
  zSVWidth_mm = 3;
  xSVStart_mm=[0 0];
  zSVStart_mm = [35 60];
  usem5 = 1;
  usem5_64 = 1; % use 64 central channels of central row
  forceRecord=0;
  largeSD=0;
  autoVSX=1;
  regionOverlapFactor3ray=1.1; % PData region overlap
  regionOverlapFactor4ray=1.2; % PData region overlap
  regionExtensionFactor=1.2; % Overall PData mag factor vs
                             % txd width
  SDLarge.zSDStart_mm = 40;
  SDLarge.zSDEnd_mm = 80;
end

if 1 | VSproductionFlag  
  productionsettings;
%  largeSDParms.SDLargeDisplayOverlay=0;
end

% override
%maxPRF_Hz = 2270;
PRIAcq_uSec = round(1e6/maxPRF_Hz);
PRIAcqInPlace_uSec = round(1e6/maxPRF_Hz);
%PRIAcq_uSec = 402;
%PRIAcqInPlace_uSec = 402;
%SDLarge.maxSpgArea_mm2 = 30*20/2;
SDLarge.zSDStart_mm = 30;
SDLarge.zSDHeight_mm = 45;
SDLarge.zSDEnd_mm = SDLarge.zSDStart_mm+SDLarge.zSDHeight_mm;
SDLarge.xSDWidth_mm = 25;

%uscmsg(execState);
%largeSDParms.largeSDSave=1;

% set window names based on version info
evParms.UI.USImageFigureTitle = [version.name ' ultrasound image window'];
evParms.UI.UCCTitle = [version.name ' ultrasound control console'];
evParms.SD.titlePrefix = [version.name ' spectrogram at '];
evParms.SD.titlePower{1} = [evParms.SD.titlePrefix 'blue' ...
                            evParms.SD.titleSuffix];
evParms.SD.titlePower{2} = [evParms.SD.titlePrefix 'green' ...
                            evParms.SD.titleSuffix];

evParms.SD.titleFlow{1} = [evParms.SD.titlePrefix ...
                    'yellow' evParms.SD.titleSuffix];
evParms.SD.titleFlow{2} = evParms.SD.titlePower{2};

evParms.UI.MImageFigureTitle = [version.name ' ultrasound M-mode window'];
evParms.UI.MImageFigureHandle = [];
evParms.UI.MImageImageHandle = [];  
evParms.UI.MImageImageGateLineHandle = [];  

if phantomFlag
  productionsettings;
  VInit = [12 12];
  % maximum voltage allowed [B/SD CFI] 
  VMax = [50 50];
  xSVWidth_mm = 5;
  zSVWidth_mm = 3;
end

if fieldMeasurementFlag 
  % settings for hydrophone measurements
  uscmsg(execState);
  productionsettings;
  disp('*** settings for hydrophone measurements with M5Sc-D');
  UTA260MUX=0;
  usem5_64=0;
  autoVSX=1;
  sendTriggers=1; % also sync at start of acq
  trig.frameNos = 2;
  trig.onB = 0;
  trig.onSD= 1;
  trig.onCFI= 0;
  trig.onBPreCFI=0;
  trig.onFirstPulseOfTypeOnly = 0;
  %%%%% SUPRATIK NOTE THIS IS ALTERNATIVE TO shift-mouse click,
  %should be same, pls check
  SDOnlyInPlace=1;
end

Trans.maxHighVoltage = max(VMax);

%largeSDParms.SDLargeCDIMode =  'freq';

% set to 1 to start in simplex high PRF mode
SDOnlyInPlace=0; 

if SDOnlyInPlace
  matFileSDOSuffix = '_SDO1';
else
  matFileSDOSuffix = '';
end

%largeSDParms.endPRI = 25;  
evParms.state.initFigsDone = 0; % ultrasound system figs not
                                % initialized yet
evParms.state.recording = [0 0]; % do not start recording IQ, per
                                 % channel
evParms.state.recordingWT = [0]; % do not start recording WT

evParms.state.seqState = seqState;
evParms.state.execState = execState;
if strcmp(evParms.state.seqState, 'BSD') 
  switch largeSDParms.SDLargeCDIMode
   case 'freq'
    evParms.state.CDIState = 3; % overlay freq
   case 'power'
    evParms.state.CDIState = 4; % overlay power
  end
else
  if strcmp(evParms.state.seqState, 'BCDI')
    % BCDI not yet supported
    uscmsg(execState);
    switch largeSDParms.SDLargeCDIMode
     case 'freq'
      evParms.state.CDIState = 1;
     case 'power'
      evParms.state.CDIState = 2;
    end
  else
    uscmsg(execState);
  end
end

SDLarge.zSDHeight_mm = SDLarge.zSDEnd_mm-SDLarge.zSDStart_mm;

% set field-of-view. Note, system will adjust these to nearest
% inclusive 128-sample boundary
% B mode is defined on a sector space for M5Sc-D and 6SD
% SD is defined on rectangle
% CFI is defined on parallelogram


FOV_B = [];
FOV_B.startDepth_mm = startDepth_mm;
FOV_B.maxDepth_mm = 85;
FOV_B.maxHalfWid_mm = 12;

FOV_SD = [];
FOV_SD.startDepth_mm = startDepth_mm;
FOV_SD.maxDepth_mm = 90;
FOV_SD.maxHalfWid_mm = 20;

FOV_CFI.startDepth_mm = SDLarge.zSDStart_mm;
FOV_CFI.maxDepth_mm = 55;
FOV_CFI.maxHalfWid_mm = 12;
  
% idea is FOV determines image size
FOV = [];
FOV.startDepth_mm = min([FOV_B.startDepth_mm ...
                    FOV_SD.startDepth_mm]);
FOV.maxDepth_mm = max([FOV_B.maxDepth_mm ...
                    FOV_SD.maxDepth_mm]);
FOV.maxHalfWid_mm = max([FOV_B.maxHalfWid_mm FOV_SD.maxHalfWid_mm]);

FOV.scanHalfAngle_deg = atand(FOV.maxHalfWid_mm/FOV.maxDepth_mm);
FOV.scanHalfAngle = FOV.scanHalfAngle_deg*pi/180;
%FOV.r1ValDop_mm = 10; % start depth CFI and SD
%FOV.r2ValDop_mm = FOV.maxDepth_mm ; % end depth CFI and SD
FOV.dopAngle = 0 * pi/180; % angle for Doppler flash transmits.
FOV.SDRegionExtensionFactor = 0; % NOT USED fraction times width of aperture

if doCFI
  % doCFISeq tries to perform CFI after B and SD. Does not work
  % well, so it is not supported
  doCFISeq=0; % cfi pri's after executing B and SD
  
  % if CFI elements are not separate subsequence, cannot have a
  % duplex mode
  if ~doCFISep
    duplexMode=0;
  end
  % use Verasonics optimized code for CDI estimation
  doCFIInternalProc=1;
else
  % do not prepare a CFI subsequence that can be switched to when needed
  doCFISeq=0;
end

if doCFI 
  if ~doCFIInternalProc
    % CFI estimation state for Matlab-based CDI code. Not supported.
    uscmsg(execState);
    evParms.state.DopStateExternalProc = 'computeCFIFreqEst';
  else
    % start in power mode
    evParms.state.DopState = 'power'; % freq
  end
end



if ~doCFISeq & ~doCFISep
  % in this case, we don't need to change the transmit or receive profile
  % because we are running B+SD sequence and there is no time to
  % switch profiles
  singleRcvProfile = 1; % keep single receive profile for B/SD/CFI
                        % (saves 2 ms/frame)
  singleTpcProfile = 1; % keep single receive profile for B/SD/CFI
                        % (saves 2 ms/frame)
else
  singleRcvProfile = 0; 
  singleTpcProfile = 0; 
end

if doCFI & ~doCFISeq & ~doCFISep
  priSkip=3; % triplex mode, not supported
else
  priSkip=2; % duplex mode, with separate B+CFI subsequence 
end

% set the Doppler scan angle for the sector probes
if (usec15 | usem5 | use6s) & doCFI
  dopScanAngle = FOV.scanHalfAngle_deg*pi/90; 
end

% settings for unsupported probes
if use312
  Trans.name = 'GEL3-12D';
  Trans.frequency = 7.8125;
  Trans.frequency = 4.4643;  
  numDopplerTx = 128;
end

if use6s
  Trans.name = 'GE6SD';
  Trans.frequency = 4.4643;
  numElementsUsed = 96;
  numDopplerTx = 96;
end

if usec15
  if usec16
    Trans.name = 'GEC1-6D';
  else
    Trans.name = 'GEC1-5D';
  end
  
  Trans.frequency = 3.4722;
  Trans.frequency = 1.9531;  
  
 uscmsg(execState);  numElementsUsed = 64; 
  numDopplerTx = 48; % approx 17mm %elements, since optic 
                       % foramen subtends about
                       % 17mm at surface of globe
end

if usec15 | use6s | use312
  uscmsg(execState);
  % for use with 256 channel Vantage system, not supported by product
  Resource.Parameters.numTransmit = 256;           % number of transmit channels.
  Resource.Parameters.numRcvChannels = 256;        % number of
                                                   % receive channels.
end

if usem5
  if usem5_64
    % settings for supported probe
    Trans.name = 'GEM5ScD_64';
    numElementsUsed = 64;
    numDopplerTx = 64; % 64 would be approx 17mm
    Resource.Parameters.numTransmit = 64; 
    Resource.Parameters.numRcvChannels = 64;
  else
    % not a formally supported configuration
    disp('*** Warning: usem5_64=0 ***');
    uscmsg(execState);
    Trans.name = 'GEM5ScD';
    %numElementsUsed = 80;
    %numDopplerTx = 80; % 64 would be approx 17mm
    numElementsUsed = 64;
    numDopplerTx = 64; % 64 would be approx 17mm
    Resource.Parameters.numTransmit = 256;
    Resource.Parameters.numRcvChannels = 256;
  end
end
    
endDepth_mm = FOV.maxDepth_mm;
maxDepth_mm = FOV.maxDepth_mm;  % maxDepth for RangeChange and RcvBuffer

% debug mode
if shallowDepth
  endDepth_mm = 45;
  maxDepth_mm = 45;
end

if doCFI
  % set focus of wide beams of CFI mode
  txFocusMm  = 4*endDepth_mm;
end

% SD and CFI demodulation frequency. It should be one of the
% ones in the Verasonics table e.g., 2.8409, 2.7174, 2.6042, 2.5000
demodFreq = Trans.frequency;

% fill evParms so we can conveniently pass the parameters to helper functions
% call figsizefn to get figure sizes
evParms.figSize = figsizefn;
evParms.ev.VInit = VInit;
% for safety, limit UI to the initialization value
evParms.ev.VMax = VMax;
evParms.ev.VMin = 1.6; % Vantage minimum voltage allowed

evParms.FOV = FOV;

if use312 | use6s |  usem5_64
  % these probes need an aperture number specified
  evParms.flag.useAperture = 1;
else
  evParms.flag.useAperture = 0;
end

evParms.flag.forceHighPRFModeForRecordingIQ= ...
    forceHighPRFModeForRecordingIQ;

% number of half cycles in excitation waveforms
evParms.TWSpec.halfCyclesB=halfCyclesB;
evParms.TWSpec.halfCyclesSD=halfCyclesSD;
evParms.TWSpec.halfCyclesCFI = halfCyclesCFI;

largeSDParms.BImGainScale = ...
    evParms.TWSpec.halfCyclesB/evParms.TWSpec.halfCyclesSD;
  
evParms.flag.separateBAcq = separateBAcq;

evParms.flag.offline=offlineSDProcFlag; % for SD offline processing, always zero
                                        % for acq scripts
evParms.flag.forceRecord = forceRecord;
evParms.ev.numSaveFramesPerBatch = numSaveFramesPerBatch;

% for saving SDLarge
iqFileSDLargePre = [vtData '/iqdata/' mfile '_' dateStr ];
evParms.file.iqFileSDLargePre = replacechar(iqFileSDLargePre, '/', '\');
evParms.file.SDLargeFramesPerFile = SDLargeFramesPerFile;
evParms.flag.startInBCFI= startInBCFI; % flag used for sdcall set

if startInBCFI
  currentModeBCFI=1;
else
  currentModeBCFI=0;
end

evParms.flag.currentModeBCFI= currentModeBCFI;
evParms.flag.audioOn = 1;
evParms.flag.instanceNoWithAudioEnabled=1;
evParms.flag.doCFINow=doCFISepAtStart;
evParms.flag.doCFI = doCFI;
evParms.flag.doCFISeq = doCFISeq;
evParms.flag.doCFISep = doCFISep;
evParms.gate.narrowbanding = narrowbanding;
evParms.gate.SDLarge = SDLarge;
evParms.gate.SDGateDisplayShape = SDGateDisplayShape;
evParms.flag.doCFISepAtStart = doCFISepAtStart;
evParms.ev.nfrms = numFrames; % no. of frames to acquire 
                              % in RcvBuffer (sets size of RF cineloop).

evParms.gate.largeGateIQ = 1; % use the larger PData for SD, as
                              % specified by {xz}SVWidthIQ_mm

evParms.flag.doWT = 0; % enable wall track via separate PDATA
evParms.flag.wtCapability = 0; % make Recon, ReconInfo for wt

% make separate bline transmit for m-mode
evParms.flag.BLineWT = 0;
evParms.flag.BLineWTCapability = 0;


if ~evParms.flag.doWT & ~evParms.flag.BLineWT
  evParms.SD.enableInstances = [1 2]; % list of enabled instances, enable sdfn calls via sdwrapfn
end

evParms.wt.xWidth_mm = 5; 
evParms.wt.zHeight_mm = 40;
evParms.wt.xDelta_wvl = 0.5;
evParms.wt.zDelta_wvl = 0.1;  

% used only in unimplemented UI control mode
if 0s
%evParms.wt.zStart_mm = 35;
%evParms.wt.xOffset_mm = 0;
end

evParms.wt.zMinHeight_mm = 10;  
evParms.wt.xMinWidth_mm = 10;    

evParms.wt.nPRI = 1; % pulses per frame

%MARK
if useCheckpoints
  disp('Before first call to SD for scheduling');
  pausede
end

evParms.ev.framePeriod_s = framePeriod_s;
evParms.ev.nPRIs = floor(2*framePeriod_s/PRIAcqInPlace_uSec*1e6)/2;
dopPRF = (1e6/PRIAcqInPlace_uSec)/2; % always the duplex value
evParms.ev.dopPRFImagingPriority = dopPRF; % if we switch back into
                                      % duplex B-SD (not SDOnlyInPlace)
                                      % later, this is used

if SDOnly | ~evParms.flag.separateBAcq
  ttna_microsec = PRIAcq_uSec;
  dopPRFUsed = dopPRF*2;
else
  ttna_microsec = PRIAcq_uSec;
  dopPRFUsed = dopPRF;  
end

evParms.ev.ttna_microsec = ttna_microsec;
evParms.ev.acqPRF = 1e6/evParms.ev.ttna_microsec;
evParms.ev.dopPRF = dopPRFUsed;
evParms.ev.priSkip = priSkip;
evParms.ev.nDopPRIs = evParms.ev.nPRIs/2;
evParms.ev.maxPRIs = maxPRIs; % maximum allocation of Receive's
evParms.ev.frameTimeActual_s = 1e-6*evParms.ev.maxPRIs*evParms.ev.ttna_microsec;
evParms.ev.PRIAcqInPlace_uSec = PRIAcqInPlace_uSec;
evParms.ev.dopPRFInPlace = 1e6/PRIAcqInPlace_uSec;


if evParms.flag.BLineWT
  % this is number of ray lines to use, currently fixed at 1
  % note: we cannot use this when switched to SD only mode
  evParms.BLineWT.tissueAtten_dBpMHz = 1.0;
  evParms.BLineWT.numChannels = 9; % channels to average around origin
  evParms.BLineWT.originTX_wvl = 0; % overwritten by gate 1 x pos
                                    % below
  evParms.BLineWT.numOrigin = length(evParms.BLineWT.originTX_wvl);
  evParms.BLineWT.numTxPerEnsemble = 20; % number of sequential
                                        % pulses at "one time point"
  evParms.BLineWT.numTxBWT = evParms.BLineWT.numOrigin*evParms.BLineWT.numTxPerEnsemble;
  evParms.BLineWT.numElTxBWT= 64; % elements to use to form focused rayline  
  evParms.BLineWT.frequency_MHz = BModeFrequency;
  TWIndBWT=4;
  evParms.ind.TWIndBWT = TWIndBWT;
  evParms.BLineWT.lambda_mm = ...
     1e-3*Resource.Parameters.speedOfSound/ ...
     evParms.BLineWT.frequency_MHz;
  lambdaTrans_mm = 1e-3*Resource.Parameters.speedOfSound/Trans.frequency;
  evParms.BLineWT.focus_mm = 65;
  evParms.BLineWT.focus_wvl = evParms.BLineWT.focus_mm/ ...
      lambdaTrans_mm; % think need to use Trans.freq value
                              % for this
  evParms.BLineWT.MMaxDepth_mm = 75;
  evParms.BLineWT.MStoreDuration_s = 60;
  evParms.numBWTTxEnsemblePerFrame = 2; % near equally spaced 
  evParms.BLineWT.frameRate_Hz = 1/evParms.ev.frameTimeActual_s;
  evParms.BLineWT.lineRate_Hz = evParms.numBWTTxEnsemblePerFrame/...
      evParms.ev.frameTimeActual_s;
  % each frame may contain more than one ensemble
  evParms.BLineWT.numFrameStore = ceil(evParms.BLineWT.MStoreDuration_s*evParms.BLineWT.frameRate_Hz);
  % B Tx's that can be replaced with BWT
  evParms.BLineWT.TxBWTStart = 5;
  evParms.BLineWT.TxBWTEnsemble = ...
      evParms.BLineWT.TxBWTStart:2:evParms.BLineWT.TxBWTStart+2*(evParms.BLineWT.numTxPerEnsemble-1);
  
  TxBWTNextFrameStart = evParms.ev.nPRIs+evParms.BLineWT.TxBWTStart;
  TXIndBWTPre = linspace(evParms.BLineWT.TxBWTStart,  ...
                         TxBWTNextFrameStart,  ...
                         evParms.numBWTTxEnsemblePerFrame+1);
 
  evParms.BLineWT.TXIndBWTPre = 2*floor(TXIndBWTPre(1:end-1)/2)+1;
  evParms.BLineWT.TxBWTInterval =  evParms.BLineWT.TXIndBWTPre(2)-...
      evParms.BLineWT.TXIndBWTPre(1);
    
  evParms.BLineWT.TXIndBWT = [];
  for q = 1:length(evParms.BLineWT.TXIndBWTPre)
    evParms.BLineWT.TXIndBWT = [evParms.BLineWT.TXIndBWT ...
                        evParms.BLineWT.TxBWTEnsemble+(q-1)* evParms.BLineWT.TxBWTInterval];
  end
  
end

% used by VSX UI sliders
minDopPRF = 1000;
maxDopPRF = dopPRF;

clear dopPRF priSkip acqPRF framePeriod_s dopPRFUsed

% samples in FFT window for SD (not including zero padding):
evParms.SD.NWindow = round(evParms.ev.dopPRF*0.014)*2; 

sgp= feval(@sdfn,'spectrogramParameters',...
           evParms.ev.framePeriod_s, evParms.SD.dopSweepTimeIndex);
evParms.SD.sweepTime_s = sgp.THist;

% set CFI widebeam parameters
if doCFI
  if use312
    % idea is to use 4 MUX apertures for L3-12D. No need for this
    % with 6SD
    dopNumRays = 4; % no. of rays in frame
  else
    dopNumRays = 3; % no. of rays in frame
  end

  if usec15 | usem5 | use6s
    % thickness of each CFI ray
    dopDtheta = dopScanAngle/dopNumRays;
  end
end

Trans.units = 'wavelengths'; 
Trans = computeTrans_nicpseq(Trans);  % computeTrans is used for known transducers.
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
lambda_m = Resource.Parameters.speedOfSound / (Trans.frequency*1e6); 
lambda_mm = lambda_m*1e3;

% setting for UTA-260MUX front panel for Vantage software version < 4.0.0
if UTA260MUX==1;
  Resource.System.UTA = '260-S'
  Trans = computeUTAMux64(Trans);
  aperture = 1;
end

% setting for UTA-260MUX front panel for Vantage software version >= 4.0.0
if UTA260MUX==2;
  Trans = computeUTAMux64(Trans);
  aperture = 1;
end

if use6s
  numElementsUsed=Trans.numelements;
end

if usec15 | usem5 | use6s
  if doCFI
%    evParms.FOV.dopR1 = FOV.r1ValDop_mm*scaleToWvl;
%    evParms.FOV.dopR2 = FOV.r2ValDop_mm*scaleToWvl;
  end
  
  if usem5
    txIndex = 1:Trans.numelements; % use only central row
    if usem5_64
      txUsedInd = txIndex;
      else
        txIndexFirst = txIndex(1:Trans.numelements/2);
        txUsedInd = findcentraln(numElementsUsed, txIndexFirst);
    end    
  end
  
  if usec15 | use6s
    txIndex = 1:Trans.numelements;
    txUsedInd = findcentraln(numElementsUsed, txIndex);
  end
end

% P contains details of the imaging geometry. It is not used by
% VSX. Once defined, it is moved into evParms
P=[];
% P index for B mode
PIndB = 1;
evParms.ind.PIndB=PIndB;
P(PIndB).startDepth = FOV_B.startDepth_mm / lambda_mm;
P(PIndB).endDepth = FOV_B.maxDepth_mm / lambda_mm;
P(PIndB).maxDepth = P(PIndB).endDepth;

if doCFI
  txFocusCFI = txFocusMm*scaleToWvl; 
  % P index for SD mode when CFI is accommodated in sequence
  PIndSD = 3;
else
  % P index for SD mode when CFI is not accommodated 
  PIndSD = 2;
end

evParms.ind.PIndSD = PIndSD;

P(PIndSD).startDepth = FOV_SD.startDepth_mm / lambda_mm;
P(PIndSD).endDepth = FOV_SD.maxDepth_mm / lambda_mm;

if focussedDoppler
  P(PIndSD).txFocus = txFocus_mm / lambda_mm;
  P(PIndSD).numTx = numDopplerTx;
end

if doCFI
  % P index for CFI mode
  
  PIndCFI=2;
  P(PIndCFI).startDepth = FOV_CFI.startDepth_mm/lambda_mm;
  P(PIndCFI).endDepth = SDLarge.zSDStart_mm/lambda_mm;
  P(PIndCFI).maxDepth = FOV_CFI.maxDepth_mm/lambda_mm;
  P(PIndCFI).dopNumRays = dopNumRays; % no. of rays in frame

  if ~doCFIInternalProc    
    % needed if Matlab code is used for CDI, not supported in product
    P(PIndCFI).DopState = evParms.state.DopStateExternalProc;
  end
else
  PIndCFI=[];
end

if doCFI
  if usec15 | usem5 | use6s
    % only single TX.aperture 
    dopTxOrgChnlPre = round([1:dopNumRays]*numElementsUsed/(dopNumRays+1));
    P(PIndCFI).dopTxOrgChnl = txUsedInd(dopTxOrgChnlPre);
    P(PIndCFI).dopNumTx = numElementsUsed;  % no. of elements in TX
                                            % aperture
                                            % each doppler beam
    % make the beams overlap, and extend reconstruction region
    % beyond the lateral extent of the transducer
    
    P(PIndCFI).regionExtensionFactor = regionExtensionFactor; 
    P(PIndCFI).regionOverlapFactor = regionOverlapFactor4ray; 
    
    if usem5_64                                            
      P(PIndCFI).dopDispEle = 64; % middle 64 elements width will
                                  % be used for
      
      P(PIndCFI).dopNumTx=Trans.numelements/2;

      if P(PIndCFI).dopNumRays==3
        % for 3 rays, beam regions do not overlap at all
        P(PIndCFI).dopNumTx = round(Trans.numelements*3/4);
        P(PIndCFI).regionOverlapFactor = regionOverlapFactor3ray;
      end
      
    else
      % not a supported case
      P(PIndCFI).dopDispEle = 128; % middle 128 elements width 
                                   % will be used for
                                   % Doppler display
    end  
  end
    
  if use312 % not a supported case
    % TX.aperture will be changed for each ray, so use central
    % element of each MUXd aperture (simpler)
    P(PIndCFI).dopTxOrgChnl = 64*ones(1,P(PIndCFI).dopNumRays); 
    % each doppler beam
    P(PIndCFI).dopDispEle = 128; % middle 128 elements width will be used for
                                 % Doppler display
    P(PIndCFI).dopNumTx = 128;  % no. of elements in TX aperture
  end
 
  if use6s % not a supported case
    uscmsg(execState);
    % only single TX.aperture 
    P(PIndCFI).dopTxOrgChnl = [1:dopNumRays]*Trans.numelements/(dopNumRays+1);
    P(PIndCFI).dopDispEle = 128; % make CFI full width
    if dopNumRays > 1
      P(PIndCFI).dopNumTx = Trans.numelements/2; % use 1/2 in 6SD,
                                                 % when multiple
                                                 % beams are used
      P(PIndCFI).dopNumTx = Trans.numelements*2/3; % use 1/2 in 6SD,
    else
      P(PIndCFI).dopNumTx = Trans.numelements;
    end
  end

  P(PIndCFI).txFocus = txFocusCFI;
  
  % put these at rate of pulses for B/SD. That means twice the
  % Doppler PRF

  if doCFISeq 
    uscmsg(execState);
    % not a supported case
    P(PIndCFI).dopPRF = evParms.ev.dopPRF*2; % Doppler PRF in Hz. (between dopPRIs tx's)
  else
    % supported case
    P(PIndCFI).dopPRF = evParms.ev.dopPRF; % Doppler PRF in Hz. (between
                                           % dopPRIs tx's)
  end

  if doCFISep
    % supported case 
    P(PIndCFI).dopStartDepth = FOV_CFI.startDepth_mm/lambda_mm;
    P(PIndCFI).dopEndDepth =   FOV_CFI.maxDepth_mm/lambda_mm;
    P(PIndCFI).dopPRIs = numCFIPRIs; 
    P(PIndCFI).dopPRF = dopPRFCFI; % Doppler PRF in Hz. 
    P(PIndCFI).pwrThres = pwrThresCFI;
    P(PIndCFI).persfreq = persfreqCFI; 
    P(PIndCFI).perspwr = perspwrCFI; 
    P(PIndCFI).cpl = cplCFI;
  end
  
  m = P(PIndCFI).dopNumRays;
end

% PData contain the formal reconstruction geometries. Used by VSX.
PDataIndB = 1;
if doCFI
  PDataIndSDStart = 4; % make room for two CFI PDatas
else
  PDataIndSDStart = 2;
end


evParms.ind.PDataIndWT = PDataIndSDStart+2+largeSD;
evParms.wt.MStoreDuration_s = 40;


PData=[];

% B-mode 
if usem5 | use6s 
  % use sector geometry for B-mode for these probes
  aperture = numElementsUsed*Trans.spacing; % aperture in wavelengths
  dapex = (aperture/4)/tand(FOV.scanHalfAngle_deg); % dist. to virt. apex, with half aperture at array.
  radius = dapex; % this just offsets us from starting region right
                  % at apex
  % Specify PData structure array.
  PData(PDataIndB).PDelta = [Trans.spacing, 0, 0.5];  % x, y and z pdeltas
  sizeRows = 10 + ceil((P(PIndB).endDepth + radius - ...
                        (radius * cosd(FOV.scanHalfAngle_deg)))/PData(PDataIndB).PDelta(3));
  sizeCols = 10 + ceil(2*(P(PIndB).endDepth + radius)*...
                       sind(FOV.scanHalfAngle_deg)/PData(PDataIndB).PDelta(1));
  PData(PDataIndB).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
  PData(PDataIndB).Origin(1,1) = (P(PIndB).endDepth+radius)*sind(-FOV.scanHalfAngle_deg) - 5; 
  PData(PDataIndB).Origin(1,2) = 0;
  PData(PDataIndB).Origin(1,3) = ceil(radius * cosd(FOV.scanHalfAngle_deg)) ...
      - radius - 5;
  PData(PDataIndB).Region = struct(...
                     'Shape',struct('Name','Sector',...
                     'Position',[0,0,-radius],...
                     'r1',radius+P(PIndB).startDepth,...
                     'r2',P(PIndB).endDepth+radius,...
                     'angle', FOV.scanHalfAngle_deg*pi/90));
  PData(PDataIndB).Region = computeRegions(PData(PDataIndB));
end

if use312
  uscmsg(execState);
  % x, y, z
  PData(PDataIndB).PDelta = [Trans.spacing, 0, 0.5];
  %z
  PData(PDataIndB).Size(1,1) = ...
      ceil((P(PIndB).endDepth-P(PIndB).startDepth)/PData.PDelta(3));
  %x
  PData(PDataIndB).Size(1,2) = ...
      ceil((Trans.numelements*Trans.spacing)/PData(PDataIndB).PDelta(1)); % cols
  %y 
  PData(PDataIndB).Size(1,3) = 1; % single image page
                                  % x,y,z of upr lft crnr w respect 
                                  % to xducer origin.
  PData(PDataIndB).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0, ...
                      P(PIndB).startDepth]; 
end

if usec15 
  uscmsg(execState);
  radius=Trans.radius;
  scanangle = numElementsUsed*Trans.spacing/radius;
  theta = -(scanangle/2); % angle to left edge from centerline

  % Specify PData structure array.
  PData(PDataIndB).PDelta = [1.0, 0, 0.5];  % x, y and z pdeltas
  sizeRows = 10 + ceil((P(PIndB).endDepth + radius - ...
                        (radius * cos(scanangle/2)))/PData(PDataIndB).PDelta(3));
  sizeCols = 10 + ceil(2*(P(PIndB).endDepth + radius)*...
                       sin(scanangle/2)/PData(PDataIndB).PDelta(1));
  PData(PDataIndB).Size = [sizeRows,sizeCols,1];     % size, origin and pdelta set region of interest.
  PData(PDataIndB).Origin(1,1) = (P(PIndB).endDepth+radius)*sin(-scanangle/2) - 5; 
  PData(PDataIndB).Origin(1,2) = 0;
  PData(PDataIndB).Origin(1,3) = ceil(radius * cos(scanangle/2)) - radius - 5;
  PData(PDataIndB).Region = struct(...
      'Shape',struct('Name','Sector',...
                     'Position',[0,0,-radius],...
                     'r1',radius+P(PIndB).startDepth,...
                     'r2',radius+P(PIndB).endDepth,...
                     'angle',scanangle));
  PData(PDataIndB).Region = computeRegions(PData(PDataIndB));
end


% - Doppler PData structure for the numSD gates
PDataIndSDEnd = PDataIndSDStart+numSD-1;
for PDataIndSD =  PDataIndSDStart:PDataIndSDEnd
  PData(PDataIndSD).PDelta = [1,0,1];
end

% TW is a structure required by VSX to define transmission waveforms
TWIndB = 1;
evParms.ind.TWIndB = TWIndB;
if doCFI
  TWIndCFI=2;
  TWIndSD=3;
else
  TWIndSD=2;
end

% coded excitation
if chirpOn
  % pulse compression
  uscmsg(execState);

  if use312
    chirp = load('ge_l3_12_d_chirp_3p5to10', 'TW');
  end
  
  if use6s
    chirp = load('ge_6sd_chirp_2p5to6p5', 'TW');
  end
  
  if usec15
    chirp = load('ge_c1_5_d_chirp2to4p8', 'TW');
  end
  
  TW(1) = chirp.TW(1); 
  TW(TWIndSD) = chirp.TW(1); 
else
  % standard pulses
  demodFreqSD = demodfreqcheckfn(Trans.frequency);
  TW(TWIndSD).type = 'parametric';
  TW(TWIndSD).Parameters = [demodFreqSD, 0.67, evParms.TWSpec.halfCyclesSD, ...
                      1];
  chirp.TW = TW;
end

% B-mode transmit waveform
TW(TWIndB).type = 'parametric';
TW(TWIndB).Parameters = [BModeFrequency,0.67, evParms.TWSpec.halfCyclesB, ...
                    1];

if evParms.flag.BLineWT
  TW(TWIndBWT).type = 'parametric';
  TW(TWIndBWT).Parameters = [evParms.BLineWT.frequency_MHz, 0.67, ...
                             evParms.TWSpec.halfCyclesB, ...
                             1];
end

evParms.ind.TWIndB = TWIndB;
evParms.ind.TWIndSD = TWIndSD;

% CFI Doppler transmit waveform, transmit frequency should be equivalent to
% supported demodFrequency
TWOffset=TWIndB+1;
if doCFI
  TWIndCFI = TWOffset;
  TW(TWIndCFI).type = 'parametric';
  evParms.ev.demodFreqCFI = demodfreqcheckfn(Trans.frequency);
  TW(TWIndCFI).Parameters = [evParms.ev.demodFreqCFI, 0.67, ...
                      evParms.TWSpec.halfCyclesCFI, 1];   
  TWOffset=TWOffset+1;
end
    
if doCFI
  PDataIndCFI1 = 2;
  PDataIndCFI2 = 3;
  if usem5 | use6s
    evParms.ind.PDataIndB=PDataIndB; 
    evParms.ind.PDataIndCFI1=PDataIndCFI1; 
    evParms.ind.PDataIndCFI2=PDataIndCFI2; 
    evParms.ind.PIndCFI=PIndCFI;
    imageWidth =2*PData(PDataIndB).Region.Shape.r2*...
                sin(PData(PDataIndB).Region.Shape.angle/2);    
    imageWidth_mm = imageWidth*lambda_mm;
    P(evParms.ind.PIndB).imageWidth_mm = imageWidth_mm;
    P(evParms.ind.PIndB).lambda_mm = lambda_mm;
    P(PIndCFI).dopAngle=FOV.dopAngle;
%    P(PIndCFI).dopR1= evParms.FOV.dopR1;
%    P(PIndCFI).dopR2= evParms.FOV.dopR2;    
    P(PIndCFI).sameShapeAsB=0; % use rectangular system for CFI vs
                               % polar for B
    P(PIndB).numRays=1; % no multiangleB support yet
    P(PIndCFI).expandCFIToEdges.flag = 1;
    P(PIndCFI).expandCFIToEdges.weight = 3/4;
    P(PIndCFI).xOffset_mm = 0;
    P(PIndCFI).xWidth_mm = inf; % allow automatic limiting    
    P(PIndCFI).zHeight_mm = (P(PIndCFI).endDepth-P(PIndCFI).startDepth)*lambda_mm;
    % P is moved into evParms here
    evParms.P = P;
    clear P
    evParms.ev.numElementsUsed = numElementsUsed;
    evParms.Trans = Trans;
    evParms.ind.TWIndCFI = 2;
    evParms.TW=TW;    
    TX=[];
  end

  if use312
    uscmsg(execState);
    TxDopOrgX = ([16 48 80 112]-63.5)*Trans.spacing;
  end
  
  if usec15    
    uscmsg(execState);
    TxDopOrgX = Trans.ElementPos(P(PIndCFI).dopTxOrgChnl,1);
  end
 
  if use6s & 0
    uscmsg(execState);
    TxDopOrgX = (P(PIndCFI).dopTxOrgChnl-(Trans.numelements-1)/2)* ...
        Trans.spacing;  
  end
end

% parameters for SD gates
evParms.gate.numSD = numSD;
for q = 1:evParms.gate.numSD
  evParms.gate.SDGateMarkerHandle{q} = ''; % stores handles
                                           % to gate UI marker objects
end

evParms.gate.PDataIndB = PDataIndB;
evParms.gate.PDataIndSDStart = PDataIndSDStart;
evParms.gate.PDataIndSDEnd = PDataIndSDEnd;
evParms.gate.xSVWidth_mm = xSVWidth_mm;
evParms.gate.zSVWidth_mm = zSVWidth_mm;
% size of IQ save regions
evParms.gate.xSVWidthIQ_mm = xSVWidthIQ_mm;
evParms.gate.zSVWidthIQ_mm = zSVWidthIQ_mm;

evParms.gate.lambda_mm = lambda_mm;
xSVStart_wvl = xSVStart_mm/lambda_mm;
zSVStart_wvl = zSVStart_mm/lambda_mm;

evParms.gate.PDataIndCFI1 = PDataIndCFI1; % using this to set rect
                                          % PData, because it's not
                                          % a sector
evParms.gate.PDataIndCFI2 = PDataIndCFI2;
    
% these are the saved gates from last run
evParms.gate.gatesFile = [vtData '/sdcgates.mat'];

if exist(evParms.gate.gatesFile, 'file')
  load(evParms.gate.gatesFile, 'gateSave');
  disp('loaded gates, overriding defaults')
  unitMm = gateSave.unitMm;
  scaleToWvl = gateSave.scaleToWvl;
  evP = evParms;
  evP.gate = gateSave;
  evParms.gate.x = evP.gate.x;
  evParms.gate.z = evP.gate.z;    
  [PData, evParms] = setpdatasdfn(PData, evParms.gate.x, ...
                                  evParms.gate.z, evParms);
  evParms.gate.xSVWidth_wvl = evP.gate.xSVWidth_wvl;
  evParms.gate.zSVWidth_wvl = evP.gate.zSVWidth_wvl;        
  
  [PData, evParms] = setpdatawtfn(PData, evParms.gate.x, ...
                                         evParms.gate.z, evParms);
  
else
  [PData, evParms] = setpdatasdfn(PData, xSVStart_wvl, ...
                                      zSVStart_wvl, evParms);
end

  
if evParms.flag.BLineWT
  evParms.BLineWT.originTX_wvl = [evParms.gate.x(1)];
  % find channels to process raw RF based on origin
  evParms = findelementsrflinefn(evParms, Trans);
   
end

if doCFI & (usem5 | use6s)
  [PData, TX, evParms] = makecfipdatatxfn(PData, TX, evParms); 
end

% set buffers using Resource structure

% - RcvBuffer(1) is for 2D, and Doppler steering angle.
Resource.RcvBuffer(1).datatype = 'int16';

if ~SDOnly
  [evParms.P] = ...
      setprisfn(evParms);
%  evParms.ev.nDopPRIsUsed=[]; % phased out
%  evParms.ev.nPRIsReduced=[]; % phased out
else
  uscmsg(execState);
  evParms.ev.nDopPRIsUsed = []; %evParms.ev.nPRIs;
  evParms.ev.nPRIsReduced = []; %evParms.ev.nPRIs;
end  

% buffer size calculation

%if evParms.flag.BLineWT
%  maxFreq_MHz = max(maxFreq_MHz, ...
%                    TW(evParms.ind.TWIndBWT).Parameters(1));
%end

% Receive buffer length detemination

numberOfSamplesPerWvl = 4; % NS200BW for B mode, at Trans.freq

FOV_B.freq_MHz = TW(evParms.ind.TWIndB).Parameters(1);
FOV_B.c_mps = Resource.Parameters.speedOfSound;
FOV_B.baseSamplesPerWvl = numberOfSamplesPerWvl;
FOV_B.baseFreq_MHz = Trans.frequency;

[samplesPerEventB, FOV_B.maxDist_wvl] = recvbufferlengthcalcfn(FOV_B);

samplesPerFrameB = evParms.ev.nPRIs/2*samplesPerEventB; % maximum number of B PRIs

FOV_SD = FOV;
FOV_SD.freq_MHz = TW(evParms.ind.TWIndSD).Parameters(1);
FOV_SD.c_mps = Resource.Parameters.speedOfSound;
FOV_SD.baseSamplesPerWvl = numberOfSamplesPerWvl;
FOV_SD.baseFreq_MHz = Trans.frequency;

[samplesPerEventSD, FOV_SD.maxDist_wvl] = recvbufferlengthcalcfn(FOV_SD);
samplesPerFrameSDDuplex = evParms.ev.nPRIs/2*samplesPerEventSD; % maximum number of B PRIs
samplesPerFrameSDOnly = evParms.ev.nPRIs*samplesPerEventSD; % maximum number of B PRIs

FOV_CFI.freq_MHz = TW(evParms.ind.TWIndCFI).Parameters(1);
FOV_CFI.c_mps = Resource.Parameters.speedOfSound;
FOV_CFI.baseSamplesPerWvl = numberOfSamplesPerWvl/2; % BS100
FOV_CFI.baseFreq_MHz = Trans.frequency;

[samplesPerEventCFI, FOV_CFI.maxDist_wvl] = recvbufferlengthcalcfn(FOV_CFI);

evParms.FOV_B = FOV_B;
evParms.FOV_SD = FOV_SD;
evParms.FOV_CFI = FOV_CFI;

BModeEventsPerFramePreCFI = 2;
samplesPerFrameCFI =  samplesPerEventB*BModeEventsPerFramePreCFI +...
    m*evParms.P(PIndCFI).dopPRIs*samplesPerEventCFI; % maximum number of B PRIs

samplesPerFrameTotalSDDuplex = samplesPerFrameB + samplesPerFrameSDDuplex;
samplesPerFrameTotalSDOnly = samplesPerFrameB + samplesPerFrameSDDuplex;

samplesPerFrameMax = max([samplesPerFrameTotalSDDuplex samplesPerFrameTotalSDOnly]);
 
if 0
minimumRcvBufferRowsPerAcq = max([rowsPerFrameB rowsPerFrameSD]);  
minLambda_m = Resource.Parameters.speedOfSound/maxFreq_MHz/1e6;
samplesPerWvlMax = ceil(numberOfSamplesPerWvl  * maxFreq_MHz/Trans.frequency);
numberOfWvlRoundTrip = (maxDistanceForRcvBufferCalc_mm/1e3)*2/minLambda_m;
minimumRcvBufferRowsPerAcq = ceil(numberOfWvlRoundTrip*samplesPerWvlMax);
end

if doCFI 
  Resource.RcvBuffer(1).rowsPerFrame = samplesPerFrameMax + samplesPerFrameCFI;
else
  uscmsg(execState);
  Resource.RcvBuffer(1).rowsPerFrame = samplesPerFrameMax;
end

Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = evParms.ev.nfrms;

% interbuffer for B-mode
interBufferIndB = 1;
evParms.ind.interBufferIndB = interBufferIndB;
Resource.InterBuffer(interBufferIndB).datatype = 'complex';
Resource.InterBuffer(interBufferIndB).numFrames = 1;  % one intermediate frame needed for 2D.

% interbuffers for SD modes
if doCFI
  interBufferIndSDStart=3;
else
  uscmsg(execState);
  interBufferIndSDStart=2;
end

evParms.flag.largeSD=largeSD;
evParms.largeSDParms = largeSDParms;
PDataIndVec = [PDataIndSDStart:PDataIndSDEnd+evParms.flag.largeSD];

if evParms.flag.largeSD
  % need to make interbuffer for large SD as large as the B-mode
  % image in case user enlarges region
  PDataIndVec(end) = evParms.ind.PDataIndB;
end

evParms.ind.interBufferIndSDLarge = interBufferIndSDStart+numSD-1+evParms.flag.largeSD;
evParms.ind.interBufferIndSDLargeOut = evParms.ind.interBufferIndSDLarge+1;
evParms.ind.interBufferIndWT = evParms.ind.interBufferIndSDLargeOut+1;

for i = 1:numSD+evParms.flag.largeSD
  % internal CDI needs complex double
  if ( (i == numSD+evParms.flag.largeSD) &  ...
       ~evParms.largeSDParms.doSDLargeInternalProc ) ...
    |  (i == numSD+evParms.flag.largeSD+1)
    uscmsg(execState);
    Resource.InterBuffer(interBufferIndSDStart+i-1).datatype = ...
        'complex single';
  else
    Resource.InterBuffer(interBufferIndSDStart+i-1).datatype = ...
        'complex';
  end
  
  Resource.InterBuffer(interBufferIndSDStart+i-1).numFrames = 1;  
  % to accommodate SDOnly and SDOnlyInPlace, this needs to be maxPRIs
  Resource.InterBuffer(interBufferIndSDStart+i-1).pagesPerFrame = ...
      evParms.ev.maxPRIs;

  Resource.InterBuffer(interBufferIndSDStart+i-1).rowsPerFrame = ...
       PData(PDataIndVec(i)).Size(1);
  Resource.InterBuffer(interBufferIndSDStart+i-1).colsPerFrame = PData(PDataIndVec(i)).Size(2);
end

if evParms.flag.wtCapability
   Resource.InterBuffer(evParms.ind.interBufferIndWT).datatype = ...
        'complex';
   Resource.InterBuffer(evParms.ind.interBufferIndWT).numFrames = ...
       1; 
   Resource.InterBuffer(evParms.ind.interBufferIndWT).pagesPerFrame ...
       = evParms.wt.nPRI;
end

if doCFI
  % interbuffer for CFI Doppler reconstructions.
  interBufferIndCFI=2;
  Resource.InterBuffer(interBufferIndCFI).datatype = 'complex';
  Resource.InterBuffer(interBufferIndCFI).numFrames = 1; % one intermediate frame needed for Doppler.
  % P.dopPRIs pages per ensemble
  Resource.InterBuffer(interBufferIndCFI).pagesPerFrame = ...
      evParms.P(PIndCFI).dopPRIs;
end

% image buffer definitions
imageBufferIndB=1;
evParms.ind.imageBufferIndB = imageBufferIndB;

% - ImageBuffer for 2D image.
Resource.ImageBuffer(imageBufferIndB).datatype = 'double';
Resource.ImageBuffer(imageBufferIndB).numFrames = evParms.ev.nfrms;

% ImageBuffer for CFI Doppler image.
imageBufferIndCFI=2;
% image buffer for CFI Doppler
Resource.ImageBuffer(imageBufferIndCFI).datatype = 'double';
Resource.ImageBuffer(imageBufferIndCFI).numFrames =  evParms.ev.nfrms;

% this image buffer can be shared between CFI and SD overlay, since
% these do not function at the same time
imageBufferIndSDLarge = imageBufferIndCFI; 

% intermediate image buffer for SD overlay prelim results
imageBufferIndSDLargeIntermed = imageBufferIndSDLarge+1;
Resource.ImageBuffer(imageBufferIndSDLargeIntermed).datatype = 'double';
Resource.ImageBuffer(imageBufferIndSDLargeIntermed).numFrames =  1; 

% - DisplayWindow for 2D image
Resource.DisplayWindow(1).Title = evParms.UI.USImageFigureTitle;
Resource.DisplayWindow(1).Type = 'Matlab'; % mod for new SDK
Resource.DisplayWindow(1).pdelta = 0.35;
ScrnSize = get(0,'ScreenSize');
%DwWidth = ceil(PData(PDataIndB).Size(2)*PData(PDataIndB).PDelta(1)/Resource.DisplayWindow(1).pdelta);
%DwHeight = ceil(PData(PDataIndB).Size(1)*PData(PDataIndB).PDelta(3)/Resource.DisplayWindow(1).pdelta);

DwWidth = ceil((2*FOV.maxHalfWid_mm/lambda_mm)/ ...
          Resource.DisplayWindow(1).pdelta);
DwHeight = ceil((FOV.maxDepth_mm/lambda_mm)/Resource.DisplayWindow(1).pdelta);

%MARK

Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
%Resource.DisplayWindow(1).ReferencePt =
%[PData(PDataIndB).Origin(1),0,PData(PDataIndB).Origin(3)];  % 2D
%imaging is in the X,Z plane
Resource.DisplayWindow(1).ReferencePt = [-DwWidth/2*Resource.DisplayWindow(1).pdelta,0,0]; 

Resource.DisplayWindow(1).numFrames = 50;
Resource.DisplayWindow(1).AxesUnits = 'mm';
%Resource.DisplayWindow(1).Type = 'Verasonics';

% - spectral Doppler transmit waveform

if use312
  uscmsg(execState);
  apertureNo = 64;
else
  apertureNo = 1;
end

if use312 | use6s
  uscmsg(execState);
  Aperture = Trans.HVMux.Aperture(:,apertureNo);
  [Elements,Dummy,Channels] = find(Aperture);
  numDelays = size(Channels,1);
end

if usec15 | usem5
  numDelays = Trans.numelements;
end

apodNotSet =  zeros(1, numDelays);
apodSet =  ones(1, numDelays);
apodSetTx = apodSet;
indRecvDoppler = 1:numDelays;

if usec15 | usem5 | use6s
  apodSet = zeros(1, numDelays);  
  apodSet(txUsedInd)=1;
  apodSetTx = apodSet;
  apodSetRx = apodSet;
  indRecvDoppler = txUsedInd;
else
  uscmsg(execState);
  apodSet =  ones(1, numDelays);
  apodSetTx = apodSet;
  apodSetRx = apodSet;
end

apodSetRxDoppler = apodNotSet;
apodSetRxDoppler(indRecvDoppler) = 1;

% Specify TX structure arrays.

if multiAngleB
  % idea would be to use multiple angle plane waves to improve
  % image quality. seems recon is too slow for inclusion with
  % SD sequence
  uscmsg(execState);
  evParms.ev.lastFrameTTNA = 0; % introduce longer delay at frame end, to see
                                % if fixes low frame rate
  numTXB = evParms.ev.nPRIs/ evParms.ev.nPRIs.priSkip;
  numRXB = numTXB;
  na = numTXB; % too slow, for recon, but reduce rinums
  
  if (na > 1)
    dtheta = (24*pi/180)/(na-1); 
    startAngle = -12*pi/180; 
  else
    dtheta = 0; 
    startAngle=0; 
  end % set dtheta to range over +/- 12 degrees.
  steerVecInd = mod((1:numTXB)-1, na)+1;
  steerVecPre = startAngle:dtheta:startAngle+(na-1)*dtheta;
  steerVec = steerVecPre(steerVecInd);
else
  if evParms.ev.lastFrameTTNA~=0
    disp('*** warning, ttna of last event will be extended');
  end
  
  numTXB = 1;
  numRXB = evParms.ev.nPRIs/evParms.ev.priSkip; 
  steerVec=BSteerAngle_deg*pi/180;
end

TXIndB = TWIndB;

% B SD CFI
if doCFI 
  numTx = numTXB + 1 + m;
else
  uscmsg(execState);
  numTx = numTXB + 1;
end

if use312 | (use6s&0)
  uscmsg(execState);
  TX = repmat(struct('waveform', 1, ...
                     'Origin', [0.0,0.0,0.0], ...
                     'focus', 0.0, ...
                     'Steer', [0.0,0.0], ...
                     'Apod', apodSetTx,...
                     'aperture', apertureNo, ...                    
                     'Delay', zeros(1,Trans.numelements)), 1, numTx);
end


% put in multiangle B-mode plane wave

if multiAngleB
  uscmsg(execState);
  % - Set event specific TX attributes.
  % recvIndBVec = TXIndB:priSkip:(TXIndB+numTXB-1)*priSkip;
  TXIndBVec = TXIndB:(TXIndB+numTXB-1);
  
  for n = 1:numTXB  % na transmit events
    TX(TXIndBVec(n)).Steer = [steerVec(n), 0.0]; %[(startAngle+(n-1)*dtheta),0.0];
    TX(TXIndBVec(n)).Delay = computeTXDelays(TX(TXIndBVec(n)));
  end
else
  % - Set event specific TX attributes.
  % recvIndBVec = TXIndB:evParms.ev.priSkip:evParms.ev.nPRIsReduced;
  TXIndBVec = TXIndB; % repmat(TXIndB, 1, numRXB);
  TX(TXIndB).Origin= [0.0,0.0,0.0];
  TX(TXIndB).focus=0;
  TX(TXIndB).Apod= apodSetTx;
  TX(TXIndB).Steer = [steerVec(1), 0.0];
  TX(TXIndB).Delay = computeTXDelays(TX(TXIndB));
end

if doCFI
  TXIndCFI = numTXB+1;
  if evParms.flag.BLineWT
    TXIndSD = TXIndCFI+m+evParms.BLineWT.numOrigin;
  else
    TXIndSD = TXIndCFI+m;
  end
else
  uscmsg(execState);
  TXIndSD = numTXB+1;
end

% -- Last TX struct needed for Doppler
if UTA260MUX
  TX(TXIndSD).aperture=1;
end

if use6s | use312
  uscmsg(execState);
  TX(TXIndSD).aperture=1;
end

TX(TXIndSD).waveform = TWIndSD;
TX(TXIndSD).Steer = [FOV.dopAngle,0.0];
[dum, indCentralEl] = findclosestinvec(Trans.ElementPos(:,1), ...
                                       xDopplerBeam_mm);  
if usem5_64
  TX(TXIndSD).Origin = [(-(numElementsUsed-1)/2 + ...
                        (indCentralEl-1)+1/2)*Trans.spacing, 0.0, 0.0];
else
  % 256 channel system
  TX(TXIndSD).Origin = [0,0,0];
end

if focussedDoppler
  TX(TXIndSD).focus = evParms.P(PIndSD).txFocus;
 
  if use312 
    uscmsg(execState);
    indCentralElL=indCentralEl-63;  % central aperture offset 63/64
                                    % from each end
    indCentralElR=indCentralEl-64;  % central aperture offset 63/64 from each end    
  else
    indCentralElL=indCentralEl;
    indCentralElR=indCentralEl;    
  end
     
  lft = indCentralElL - floor(evParms.P(PIndSD).numTx/2)+1;
  if lft < 1, lft = 1; end;
  rt = indCentralElR + floor(evParms.P(PIndSD).numTx/2);
  if rt > Trans.numelements, rt = Trans.numelements; end;

  TX(TXIndSD).Apod = apodNotSet;
  TX(TXIndSD).Apod(lft:rt) = 1.0;
else 
  TX(TXIndSD).Apod = apodSet;
  TX(TXIndSD).focus = 0; % plane wave
  
end

TX(TXIndSD).Delay = computeTXDelays(TX(TXIndSD));  






% from widebeamDoppler:
% -- P.dopNumRays TX structs needed for Doppler
%winNum = round(192*Trans.numelements/128);

winNum = 128;
W = hannfn(winNum)'; % 128 elememts for hann window

if use312
  uscmsg(execState);
  aper = [17 49 81 113]; % center 128 elements will be used for
                         % doppler
else
  aper = 1; 
end

if useCheckpoints
  disp('TX complete');
  pausede
end

% note, TPCIndSD is not used, because we do not set a new TPC for
% SD, since we do not want to introduce delay by reseting TPC. but
% can we afford to?

evParms.ind.TPCIndB = 1;
evParms.ind.TPCIndCFI = 2;
evParms.ind.TPCIndSD = 1;

% Specify TPC structures. (adapted from widebeamDoppler)
TPC(evParms.ind.TPCIndB).name = '2D';
TPC(evParms.ind.TPCIndB).maxHighVoltage = Trans.maxHighVoltage;

TPC(evParms.ind.TPCIndCFI).name = 'CFI';
TPC(evParms.ind.TPCIndCFI).maxHighVoltage = Trans.maxHighVoltage;

if evParms.ind.TPCIndB == evParms.ind.TPCIndSD
  TPC(evParms.ind.TPCIndSD).name = 'B + Spectral Doppler';
else
  uscmsg(execState);
  TPC(evParms.ind.TPCIndSD).name = 'Spectral Doppler';
  TPC(evParms.ind.TPCIndSD).maxHighVoltage = Trans.maxHighVoltage;
end

% set receive profiles
% note that current sequences use the high gain SD profile also for
% B-mode

RcvProfile = [];

evParms.ind.RcvProfileIndB = 1;
evParms.ind.RcvProfileIndCFI = 2;
evParms.ind.RcvProfileIndSD = 3;

RcvProfile(evParms.ind.RcvProfileIndB).LnaGain = 18; % Profile used for B-mode imaging

if evParms.flag.BLineWT
  RcvProfile(evParms.ind.RcvProfileIndB).LnaGain = 24;
end


% Profile used for Doppler. Force high-Z state for best Doppler sensitivity
RcvProfile(evParms.ind.RcvProfileIndCFI).LnaGain = 24;
RcvProfile(evParms.ind.RcvProfileIndCFI).LnaZinSel = 31; 

RcvProfile(evParms.ind.RcvProfileIndSD).LnaGain = 24; % Profile used for Doppler
RcvProfile(evParms.ind.RcvProfileIndSD).LnaZinSel = 31; 

TGCIndB=1;
if doCFI
  TGCIndCFI=2;
  TGCIndSD=3;
  % CFI Doppler TGC waveform
  TGC(TGCIndCFI).CntrlPts = [611.3920 764.2400 978.2272 1023 1023 ...
                      1023 1023 1023]; %[400,500,640,710,770,830,890,950];
  TGC(TGCIndCFI).rangeMax = evParms.P(PIndB).endDepth;
  TGC(TGCIndCFI).Waveform = computeTGCWaveform(TGC(TGCIndCFI));
else
  TGCIndSD=2;
  TGCIndCFI=[];  
end

evParms.ind.TGCIndB = TGCIndB;
evParms.ind.TGCIndSD = TGCIndSD;
evParms.ind.TGCIndCFI = TGCIndCFI;

% Specify TGC Waveform structures.
evParms.default.TGC(TGCIndB).CntrlPts = 1023*ones(1,8);
evParms.default.TGC(TGCIndSD).CntrlPts = 1023*ones(1,8);
evParms.default.TGC(TGCIndCFI).CntrlPts = 1023*ones(1,8);
evParms.default.TGC(TGCIndB).rangeMax = evParms.P(PIndB).endDepth;
evParms.default.TGC(TGCIndSD).rangeMax = evParms.P(PIndB).endDepth;
evParms.default.TGC(TGCIndCFI).rangeMax = evParms.P(PIndB).endDepth;

TGC(TGCIndB).CntrlPts = evParms.default.TGC(TGCIndB).CntrlPts;
TGC(TGCIndB).rangeMax = evParms.default.TGC(TGCIndB).rangeMax;
TGC(TGCIndB).Waveform = computeTGCWaveform(TGC(TGCIndB));
% Spectral Doppler TGC waveform
TGC(TGCIndSD).CntrlPts = evParms.default.TGC(TGCIndSD).CntrlPts;
TGC(TGCIndSD).rangeMax = evParms.default.TGC(TGCIndSD).rangeMax;
TGC(TGCIndSD).Waveform = computeTGCWaveform(TGC(TGCIndSD));

evParms.TGC=TGC;

% Specify Receive structure arrays -
%   Define enough Receives to handle the maxPRI case.  For lower PRIs, some of the Receives will go
%   unused in the Event list.
%   endDepth - add additional acquisition depth to account for some channels
%              having longer path lengths.
%   InputFilter - The same coefficients are used for all channels. The
%              coefficients below give a broad bandwidth bandpass filter.

% from widebeamDoppler

evParms.rcvParms=[];
wl4sPer128 = 128/(4*2);  % wavelengths in a 128 sample block for 4
                         % smpls per wave round trip.
evParms.rcvParms.wl4sPer128 = wl4sPer128;  

if doCFI
  % maximum lateral extent of CFI region in wavelengths
  maxLateralCFI = evParms.P(PDataIndCFI1).maxLateralExtentCFI;
  if maxLateralCFI*lambda_mm > FOV.maxHalfWid_mm
    error(['CFI lateral extent larger than lateral FOV. Something is ' ...
          'wrong']);
  end

  maxAcqLengthCFI =  sqrt(evParms.P(PIndCFI).dopEndDepth^2 + (maxLateralCFI*2)^2) - ...
                          evParms.P(PIndCFI).dopStartDepth;
  
  samplesPerWaveDop = 2*TW(TWIndCFI).Parameters(1)/Trans.frequency;

  % BPF copied from widebeamDoppler
  % 30% BW
  BPF = ...
      [+0.00021 +0.00000 -0.00131 +0.00000 +0.00443 +0.00000 -0.01099 ...
       +0.00000 +0.02228 +0.00000 -0.03876 +0.00000 +0.05957 +0.00000 ...
       -0.08209 +0.00000 +0.10257 +0.00000 -0.11697 +0.00000 +0.12213];
  wl2sPer128 = 128/(2*samplesPerWaveDop);  % wavelengths in a 128
                                           % sample block for 2 smpls
  
  evParms.rcvParms.wl2sPer128 = wl2sPer128;
  evParms.rcvParms.maxAcqLengthCFI = maxAcqLengthCFI;
end

indRecv = find(apodSetRx);
indRecvCFI = indRecv;
indRecvDoppler = find(apodSetRxDoppler);

% Set event-specific Receive attributes.

evParms.rcvParms.PIndB = PIndB;
evParms.rcvParms.PIndSD = PIndSD;
evParms.rcvParms.PIndCFI = PIndCFI;
evParms.rcvParms.sampleMode.B = 'NS200BW';
evParms.rcvParms.sampleMode.SD = 'NS200BW';
%evParms.rcvParms.sampleMode.SD = 'BS100BW';


if exist('aperture', 'var') & ~isempty(aperture)
  evParms.rcvParms.aperture = aperture;
  evParms.rcvParms.apertureNo = apertureNo;
end

% this is badly named. it is Receive.endDepth, and 
% sets end time of last receive. It is not roundtrip.
%evParms.rcvParms.maxAcqLength2D = numberOfWvlRoundTrip/2;
evParms.rcvParms.apodNotSet = apodNotSet;
evParms.rcvParms.indRecv = indRecv;
evParms.rcvParms.indRecvDoppler = indRecvDoppler;
evParms.rcvParms.indRecvCFI = indRecvCFI;
evParms.rcvParms.BIndVal = 1; % indicates B-mode
evParms.rcvParms.SDIndVal = 2; % indicates spectral Doppler mode, jsm
                       % naming only
evParms.rcvParms.CFIIndVal = 3; % indicates CFI mode
evParms.rcvParms.TGCIndSD = TGCIndSD;
evParms.rcvParms.TGCIndB = TGCIndB;
evParms.rcvParms.TGCIndCFI = TGCIndCFI;

if doCFI
  evParms.rcvParms.BPF = BPF;
  evParms.rcvParms.TW = TW;
  evParms.rcvParms.TWIndCFI = TWIndCFI;
end

% start of each CFI ensemble
rcvCFIStartInd = [];

evParms.ev.priSkipRcv = 2; % make enough receives for SDOnly mode
evParms.state.SDOnly=SDOnly;
evParms.flag.multiAngleB=multiAngleB;
evParms.ev.numTXB = numTXB;
evParms.flag.use312 = use312;
evParms.flag.use6s = use6s;
evParms.flag.usec15 = usec15;
evParms.flag.usem5= usem5;
evParms.ind.PDataIndCFI1 = PDataIndCFI1;
evParms.ind.PDataIndCFI2 = PDataIndCFI2;
evParms.ind.interBufferIndCFI = interBufferIndCFI;
evParms.ind.PDataIndSDStart = PDataIndSDStart;
evParms.ind.interBufferIndSDStart = interBufferIndSDStart;
evParms.ev.TXIndB = TXIndB;
evParms.ev.TXIndSD = TXIndSD;

clear PDataIndCFI1 PDataIndCFI2

if doCFI
  evParms.ev.TXIndCFI = TXIndCFI;
  evParms.flag.doCFISeq = doCFISeq;
end

evParms.ev.preFrameHWDelay_s = preFrameHWDelay_s;
evParms.flag.sendTriggers = sendTriggers;
evParms.trig=trig;
evParms.flag.syncFrameStart = syncFrameStart;
evParms.flag.singleRcvProfile = singleRcvProfile;
evParms.flag.singleTpcProfile = singleTpcProfile;
evParms.ev.TXIndBVec = TXIndBVec;

% - For 2D, we need P.numRays ReconInfo structures for P.numRays
% steering angles
evParms.flag.doCFIInternalProc = doCFIInternalProc;
evParms.ind.imageBufferIndCFI =  imageBufferIndCFI; 
evParms.ind.imageBufferIndSDLarge = imageBufferIndSDLarge;
evParms.ind.imageBufferIndSDLargeIntermed = ...
    imageBufferIndSDLargeIntermed;
evParms.ev.TXIndCFI = TXIndCFI;
evParms.ev.TXOffset =  TXIndCFI-1;

Control = [];
SeqControl = [];
Receive = [];
seq = [];

% development code for sequence on-the-fly reprogramming
if 0
  uscmsg(execState);
  evParms.state.SDOnly=1;
  dopPRFNom = 4500; % this realized ttna of 222us between SD acq, which
                    % is what is needed
  dopPRFNom = 2500

  if evParms.state.SDOnly
    ttna_microsec = round(1/dopPRFNom/1e-6);
    nPRIs = round(evParms.ev.framePeriod_s/ttna_microsec*1e6);
  else
    ttna_microsec = round(1/dopPRFNom/2e-6);
    nPRIs = round(1/2*evParms.ev.framePeriod_s/ttna_microsec*1e6)*2; % make it even
  end

  if nPRIs > evParms.ev.maxPRIs
    disp('nPRIs > max. setting to max');
    nPRIs = evParms.ev.maxPRIs
    ttna_microsec = round(1e6*framePeriod_s/nPRIs)
  end

  evParms.ev.ttna_microsec = ttna_microsec;
  evParms.ev.nPRIs = nPRIs;

  if evParms.state.SDOnly % remove interleaving for SD only mode
    evParms.ev.dopPRF=1/(evParms.ev.ttna_microsec*1e-6);
    nPRIs = round(evParms.ev.framePeriod_s/ttna_microsec*1e6);
    evParms.ev.priSkip = 1;
    evParms.ev.nDopPRIs = nPRIs;
    evParms.ev.nDopPRIsUsed = nPRIs;  
  else
    evParms.ev.dopPRF=1/(evParms.ev.ttna_microsec*2e-6);
    evParms.ev.nDopPRIs = evParms.ev.nPRIs/2;  
    evParms.ev.nDopPRIsUsed = evParms.ev.nPRIs/2;    
    evParms.ev.priSkip = 2;
  end
end % development code

% development flag to suppress running of ROIPlot at start
evParms.flag.runROIPlotForSDOnly=0;
% set to 1 to start with SDOnlyInPlace. 0 is only supported option
evParms.state.SDOnlyInPlace = SDOnlyInPlace;

% 
seqContainer = makenicpseqfn(evParms, Resource, SeqControl, Receive, ...
                             seq);

evParms = seqContainer.evParms;
recon = seqContainer.recon;
ri = seqContainer.ri;
Receive = seqContainer.Receive;
Recon = seqContainer.Recon;
ReconInfo = seqContainer.ReconInfo;
Resource = seqContainer.Resource;
Process = seqContainer.Process;
Event = seqContainer.Event;
SeqControl = seqContainer.SeqControl;

% development code for testing on-the-fly changes to the sequence
if 0
  uscmsg(execState);
  disp('*****')
  evParms.state.SDOnlyInPlace = 0;
%  evParms.flag.largeSD=0;
  if evParms.state.SDOnlyInPlace
    UIValue = evParms.ev.dopPRFInPlace;
  else
    UIValue = evParms.ev.acqPRF/2;
  end
  
  changeprffn(UIValue);
end

% development code 
if 0 
  uscmsg(execState);
  
  if 0 % suppressing ultrasound display window
    %Resource = rmfield(Resource, 'DisplayWindow');
    %Process=Process(1:2);
  end
  
  if 0 % for testing for possible memory leak in VSX
    TX(1)=TX(5);
    TX=TX(1);
    ReconOrig=Recon;
    Recon(1)=ReconOrig(3);
    Recon(2)=ReconOrig(4);
    Recon= Recon(1:2);
    ReconInfo=ReconInfo(1:450);
  end

  if 0 % for testing for possible memory leak in VSX
    TXOrig=TX;
    TX(1)=TXOrig(5);
    TX(5)=TXOrig(1);
    ReconOrig=Recon;
    Recon(1)=ReconOrig(3);
    Recon(2)=ReconOrig(4);
    Recon(3)= ReconOrig(1);
    Recon(4)= ReconOrig(2);
  end

  if 0 % for testing for possible memory leak in VSX
    ReconOrig=Recon;
    Recon(1)=ReconOrig(3);
    Recon(2)=ReconOrig(4);
    Recon(3)= ReconOrig(1);
    Recon(4)= ReconOrig(2);
  end
  
  if 0 
    %ReconInfo=ReconInfo(1:450);
    %Recon(1)=Recon(3);
    %Recon(2)=Recon(4);
    %Recon= Recon(1:2);
  end

  
  % maximum number of events ever
  %evParms.ev.numEventMax = length(Event);
  
  if 0
    evParms.state.SDOnlyInPlace = 1;
    evParms.ev.nDopPRIs=evParms.ev.nPRIs;
    evParms.ev.dopPRF=evParms.ev.acqPRF;
    
    %keyboard
    UIValue = evParms.ev.dopPRF;
    changeprffn(UIValue);
  end
  
end % development code

% *** UI structure definition for VSX controls and mouse/keyboard shortcuts
% User specified UI Control Elements
% - Sensitivity Cutoff (overall). Adjusts all Recon elements in callback
if isfield(recon, 'reconIndB')
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',...
                  [0,1.0,Recon(recon.reconIndB).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%-UI#1Callback');
end

% - Range Change. 
MinMaxVal = [64,300,evParms.P(evParms.ind.PIndB).endDepth]; % default unit is wavelength
AxesUnit = 'wls';
if isfield(Resource, 'DisplayWindow')
  if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        AxesUnit = 'mm';
        MinMaxVal = MinMaxVal * (Resource.Parameters.speedOfSound/1000/Trans.frequency);
    end
  end
end

UI(2).Control = {'UserB6','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',MinMaxVal,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%-UI#2Callback');

% - PRF control
UI(3).Control = {'UserB4','Style','VsSlider','Label','Doppler PRF',...
    'SliderMinMaxVal',[minDopPRF,maxDopPRF*2,evParms.ev.dopPRF],...
    'SliderStep',[100/ (maxDopPRF-minDopPRF),200/ (maxDopPRF-minDopPRF)],'ValueFormat','%4.0f Hz'};
UI(3).Callback = text2cell('%-UI#3Callback');

% - Wall Filter control
UI(4).Control = {'UserB3','Style','VsSlider','Label','Doppler WF',...
                  'SliderMinMaxVal',[8,0.1*maxDopPRF,evParms.SD.dopWFCutoffNorm*evParms.ev.dopPRF],...
                  'SliderStep',[0.02,0.2],'ValueFormat','%4.0f Hz'};
UI(4).Callback = text2cell('%-UI#4Callback');


% - Sweep Time control
UI(5).Control = {'UserB2','Style','VsSlider','Label','Sweep Time',...
                  'SliderMinMaxVal',[min(sgp.THistAllow),max(sgp.THistAllow), evParms.SD.sweepTime_s],...
                  'SliderStep',[0.25,0.25],'ValueFormat','%1.0f s'};
UI(5).Callback = text2cell('%-UI#5Callback');

% - Baseline control
UI(6).Control = {'UserB1','Style','VsSlider','Label','Baseline',...
                  'SliderMinMaxVal',[-0.5,0.5,evParms.SD.baseLineShiftNorm],...
                  'SliderStep',[0.05,0.1],'ValueFormat','%1.2f'};
UI(6).Callback = text2cell('%-UI#6Callback');

% - Window ButtonDown callback for setting sample volume position.
enableKeyEvents=0;
if enableKeyEvents
  % allow keyboard events to be used.
  uscmsg(execState);
  UI(7).Statement = ['set(Resource.DisplayWindow(1).figureHandle,' ...
                   '''WindowButtonDownFcn'',@nicpseqwt_wbdCallback' ...
                   ', ''WindowKeyPressFcn'', @wkbdCallback);'];
else  
  UI(7).Statement = @()setCallback(Resource.DisplayWindow(1).Type);
  %UI(7).Statement = ['set(Resource.DisplayWindow(1).figureHandle,' ...
  %                   '''WindowButtonDownFcn'',@wbdCallback)'];
end
UI(7).Callback = @nicpseqwt_wdbCallback; %(hObject,eventdata);

uiInd=8;

if 0 % narrowbanding control, development only 
    uscmsg(execState);
    UI(8).Control = {'UserB5','Style','VsButtonGroup','Title','narrowbanding',...
                     'NumButtons', 2, 'Label', {'Off', 'On'}};
    UI(8).Callback = text2cell('%-UI#8Callback');
    uiInd=uiInd+1;
end

if 0 % enable key events
  uscmsg(execState);
  UI(uiInd).Statement = ['set(Resource.DisplayWindow(1).figureHandle' ...
                      ', ''WindowKeyPressFcn'', @wkbdCallback);'];
  UI(uiInd).Callback = text2cell('%-UI#9Callback');
  uiInd=uiInd+1;
end

% - CFI PRF adjustment
if doCFI
  PRFmin = 500; PRFmax = 4500; stepDiff = PRFmax-PRFmin;
  slider.prfcfi.pos = 'UserB5';
  UI(uiInd).Control = {slider.prfcfi.pos,'Style','VsSlider','Label',...
                      'CFI PRF (Hz)','SliderMinMaxVal',...
                      [PRFmin,PRFmax,evParms.P.dopPRF],...
                      'SliderStep',[50/stepDiff,500/stepDiff],...
                      'ValueFormat','%3.0f'};
  UI(uiInd).Callback = text2cell('%PRFChange');
  evParms.ind.UIIndCFIPRF = uiInd;
  
  % - Steer Angle adjustment
  % no space on VSX display to show this, overwritten by MaxPower
  uiInd=uiInd+1;
  slider.cfiangle.pos = 'UserC7';
  UI(uiInd).Control = {slider.cfiangle.pos,'Style','VsSlider','Label',...
                       'Doppler Angle','SliderMinMaxVal',...
                      [-20,20,round(evParms.P(PIndCFI).dopAngle*180/pi)],...
                      'SliderStep',[1/40,5/40],'ValueFormat','%3.0f'};
  UI(uiInd).Callback = text2cell('%SteerAngle');

  % - CFI Doppler Mode Button Group: allows selection of velocity
  % or power. Note, will change the color map for SD overlay
  % TODO: disable when in overlay mode

  if ~Mcr_GuiHide 
      button.dopplermode.pos = 'UserC4';
      UI(uiInd).Control = {button.dopplermode.pos,'Style','VsButtonGroup','Title','Doppler Mode',...
                          'NumButtons',2,'Labels',{'Velocity','Power'}};
      UI(uiInd).Callback = text2cell('%DopplerModeCallback');
      button.uiInd=uiInd; % use this to prevent rendering of this box
      uiInd=uiInd+1;
  else
      button=[];
  end
  
  
  slider.dopplerpwrthresh.pos = 'UserC3';

  % - Doppler Power Threshold Slider
  if 0 % version for control of CDI
      uscmsg(execState);
      UI(uiInd).Control = {slider.dopplerpwrthresh.pos,'Style','VsSlider',...
                          'Label','DopPwrThres','SliderMinMaxVal',[0.0,1.0,evParms.P(PIndCFI).pwrThres],...
                          'SliderStep',[0.02,0.1],'ValueFormat', ...
                          '%3.2f'};
  end
  
  % version for control of SD overlay
  UI(uiInd).Control = {slider.dopplerpwrthresh.pos,'Style','VsSlider',...
                      'Label','DopPwrThres','SliderMinMaxVal', ...
                      [0.0,0.2, evParms.largeSDParms.pwrThres],...
                      'SliderStep',[0.0025,0.005],'ValueFormat', ...
                      '%3.4f'};

  UI(uiInd).Callback = setprocparmfn('pwrThreshold', ...
                                     'procIndSDLarge', 'evParms.largeSDParms.pwrThres');
  
  uiInd=uiInd+1;
  
  slider.dopplerMaxPower.pos = 'UserC8';
  
  % - Doppler maxPower Slider
  UI(uiInd).Control = {slider.dopplerMaxPower.pos,'Style','VsSlider',...
                      'Label','MaxPower','SliderMinMaxVal', ...
                      [0,100, evParms.largeSDParms.maxPower],...
                      'SliderStep',[0.01,0.1],'ValueFormat', ...
                      '%3.0f'};

  UI(uiInd).Callback = setprocparmfn('maxPower', ...
                                     'procIndSDLarge', 'evParms.largeSDParms.maxPower');
   
  uiInd=uiInd+1;
   
  slider.colorprioritythresh.pos = 'UserC2';
  % - Color Priority Threshold Slider
    
  switch evParms.state.CDIState
   case 1
    uscmsg(execState);
    updateParm = 'evParms.P(PIndCFI).cpl';
    procIndStr = 'procIndCFIIm';
   case 2
    uscmsg(execState);
    updateParm = 'evParms.P(PIndCFI).cpl';
    procIndStr = 'procIndCFIIm';
   case 3
    updateParm = 'evParms.largeSDParms.cpl';
    procIndStr = 'procIndSDLargeIm';
   case 4
    updateParm = 'evParms.largeSDParms.cpl';
    procIndStr = 'procIndSDLargeIm';
  end
    
  cplParm = eval(updateParm);
  
  UI(uiInd).Control = {slider.colorprioritythresh.pos,'Style','VsSlider',...
                      'Label','Color Priority',...
                      'SliderMinMaxVal',[0,255,cplParm],...
                      'SliderStep',[1/255,0.1],'ValueFormat', ...
                      '%3.0f'};
  UI(uiInd).Callback = setprocparmfn('threshold', procIndStr, ...
                                                  updateParm);
  uiInd=uiInd+1;

  % - Color Persistence Slider
  slider.perspwr.pos = 'UserC1';
  slider.perspwr.tag = [slider.perspwr.pos 'Slider'];

  switch evParms.state.CDIState
   case 1 % conventional CDI
    uscmsg(execState);
    updateParm = 'evParms.P(PIndCFI).persfreq';
    procIndStr = 'procIndCFIIm';
   case 2 % conventional CDI
    uscmsg(execState);
    updateParm = 'evParms.P(PIndCFI).perspwr';
    procIndStr = 'procIndCFIIm';
   case 3 % SD overlay CDI, flow mode
    uscmsg(execState);
    updateParm = 'evParms.largeSDParms.persfreq';
    procIndStr = 'procIndSDLargeIm';
   case 4 % SD overlay CDI, power mode
    updateParm = 'evParms.largeSDParms.perspwr';
    procIndStr = 'procIndSDLargeIm';
  end
    
  persParm = eval(updateParm);
  % - Color Persistence Slider
  UI(uiInd).Control = {slider.perspwr.pos,'Style','VsSlider','Label','Color Persistence',...
                        'SliderMinMaxVal',[0,100, persParm],...
                        'SliderStep',[1/100,0.1],'ValueFormat','%3.0f'};
 
  UI(uiInd).Callback = setprocparmfn('persistLevel', procIndStr, ...
                                     updateParm);
  
  uiInd=uiInd+1;
end % doCFI

% change the button initial position to Power
if doCFI & evParms.state.CDIState == 2
  uscmsg(execState);
  UI(uiInd).Statement = ...
      ['set(findobj(''Tag'',''UserC4RadioButton2''), ''Value'','...
                      ' 1);'];
  uiInd=uiInd+1;
end

% set initial voltages
UI = inithvfn(UI, uiInd, VInit);
uiInd=uiInd+1;

%if doCFI
%  % - External function for ROIplot
%  EF(1).Function = text2cell('%-EF#1');
%end

matOutPath = veraPath;

% Save all the structures to a .mat file.
if 0 
  uscmsg(execState);
  matFileOut = [matOutPath mfile '_' dateStr '.mat'];
  save(matFileOut);
  lslrt(matFileOut);
end


%matFileOut = [mfile '.mat'];

matFileOut = [matOutPath '/' mfile matFileSuffixSim matFileSDOSuffix ...
'.mat'];

%matFileOut = [matOutPath '/nicpseq' matFileSuffixSim matFileSDOSuffix ...
%'.mat'];

evParms.UI.slider=slider;
evParms.UI.button=button;

clear vsxGUIh updateh slider msg hvSldr hvValue fname  ...
    TXEventCheckh Control eventdata SDop

% trick to run VSX automatically

filename= matFileOut;

if 1
save(matFileOut, 'Resource', 'PData', 'Event', 'Recon', 'ReconInfo', ...
     'Receive', 'TX', 'TW', 'Process', 'SeqControl', 'Trans', 'TGC', ...
     'RcvProfile', 'evParms', 'UI', 'filename');
end
%save(matFileOut);
lslrt(matFileOut);

if Resource.Parameters.simulateMode == 2 
 disp('Press any key to load simulation data');
 pausede
 disp('Loading simulation data. Can take a while ...')
  rcvDataFile = [vtData '/rf/' ...
                 'RFdata_20-January-2019_17-05-02.mat'];
  load(rcvDataFile, 'RcvData');
  disp('Loaded simulation data')
end

  
if autoVSX
  VSX
end

return
% end of nicpseq.m base section

% **** Callback routines to be encoded by text2cell function. ****

%-UI#1Callback - Sensitivity cutoff change
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
%-UI#1Callback

%-UI#2Callback - Range change
disp('Range change not yet supported');
% To do:
% 1. For CFI, reduce range if UIValue is less than CFI limit, but max
% it out at CFI limit if UIValue is larger
% 2. Depths need to be adjusted for B and SD together.
% 3. Mostly reuse existing code to do this.
return
simMode = evalin('base','Resource.Parameters.simulateMode');
% No range change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','max([evParms.P.endDepth])'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','evParms.P');
P(PIndB).endDepth = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits') && ~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        P(PIndB).endDepth = UIValue*scaleToWvl;    
    end
end
assignin('base','P',P);

evalin('base','PData(PDataIndB).Size(1) = ceil((P(PIndB).endDepth-P(PIndB).startDepth)/PData(PDataIndB).PDelta(3));');
evalin('base','PData(PDataIndB).Region = computeRegions(PData(PDataIndB));');
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(PDataIndB).Size(1)*PData(PDataIndB).PDelta(3)/Resource.DisplayWindow(1).pdelta);');
Receive = evalin('base', 'Receive');

maxAcqLength = ceil(sqrt(P(PIndB).endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC(TGCIndB).rangeMax = P(PIndB).endDepth;');
evalin('base','TGC(TGCIndB).Waveform = computeTGCWaveform(TGC(TGCIndB));');
evalin('base','if VDAS==1, Result = loadTgcWaveform(1); end');

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','Receive','Recon','DisplayWindow'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%-UI#2Callback

%-UI#3Callback - Doppler PRF for SD
disp(['PRF change not yet supported. Switch to SDOnlyInPlace to ' ...
      'increase Doppler PRF to 4.5 kHz']);
% To do:
% resolve getStruct: Event.rcv value greater than no. of Receive
% structures. Probably it is best to define maximum number of
% Receive structures. However, may still have issues with number of
% events changing.
return
changeprffn(UIValue);
%-UI#3Callback

%-UI#4Callback - Doppler Wall Filter
evParms = evalin('base','evParms'); 
dopPRF = evParms.ev.dopPRF;
evParms.SD.dopWFCutoffNorm = UIValue/dopPRF;
assignin('base', 'evParms', evParms);
return
%-UI#4Callback

%-UI#5Callback - Sweep Time
sweepTime = round(UIValue);
evParms = evalin('base','evParms'); 
TFrame=evParms.ev.framePeriod_s;
sgp=feval(@sdfn,'spectrogramParameters',TFrame,1);
thisSweepTimeAllowed = ismember(sweepTime, sgp.THistAllow);
if thisSweepTimeAllowed
  [~,swInd] = min(abs(sgp.THistAllow-sweepTime));
  sweepTimeAllowed = sgp.THistAllow(swInd);  
  dopSweepTimeIndex = swInd;
else
  currentSweepTimeIndex = evParms.SD.dopSweepTimeIndex;
  currentSweepTime_s =   sgp.THistAllow(currentSweepTimeIndex);
  if sweepTime < currentSweepTime_s
    ind = find(sgp.THistAllow <= sweepTime, 1,  'last');
    dopSweepTimeIndex = ind;
  else
    ind = find(sgp.THistAllow >= sweepTime, 1, 'first');
    dopSweepTimeIndex = ind;
  end
  if isempty(ind)
    error('Error finding a valid sweep time!');
  end
end
newSweepTime_s = sgp.THistAllow(dopSweepTimeIndex);
evParms.SD.dopSweepTimeIndex = dopSweepTimeIndex;
evParms.SD.sweepTime_s = newSweepTime_s;
assignin('base', 'evParms', evParms);
set(hObject,'Value', newSweepTime_s);
objEdit = findobj('tag','UserB2Edit');
set(objEdit, 'string', [num2str(newSweepTime_s, '%1.0f') ' s']);
return
%-UI#5Callback

%-UI#6Callback - Baseline
baseLineShiftNorm = UIValue;
baseLineShiftNorm = round(baseLineShiftNorm*20)/20;
set(hObject,'Value',baseLineShiftNorm);
evalin('base', ['evParms.SD.baseLineShiftNorm = ' ...
                num2str(baseLineShiftNorm) ';']);
return
%-UI#6Callback

%%% mouse buttondown callback

%-UI#8Callback - narrowbanding off
%narrowbandingUICallback(hObject, eventdata);
  S = get(eventdata.NewValue,'Tag');
  UIState = str2double(S(18));
  assignin('base', 'narrowbanding', UIState-1);
  evalin('base', 'narrowbanding');
return
%-UI#8Callback

%-UI#9Callback - CFI switch on B-mode
wkbdCallback(hObject,eventdata)
evParms = evalin('base', 'evParms');
% capture keyboard presses
keyEvent = get(hObject,'CurrentCharacter')
return
%-UI#9Callback 

%DopplerModeCallback - Doppler mode change
evParms = evalin('base', 'evParms');
P = evParms.P;
PIndCFI = evParms.rcvOut.PIndCFI;
slider = evParms.UI.slider;

Control = repmat(struct('Command','set&Run','Parameters',[]),1,4);
Process = evalin('base','Process');
Resource = evalin('base','Resource');
hDisplay = Resource.DisplayWindow(1).figureHandle;
currentMap = get(hDisplay,'Colormap');
switch UIState
    case 1  % Velocity mode
        newMap = grayscaleCFImap;
        newMap(1:128,:) = currentMap(1:128,:);
        P(PIndCFI).perspwr = get(findobj('Tag',slider.perspwr.tag),'Value');
        Control(1).Parameters = {'Process',evParms.proc.procIndCFI,'method','computeCFIFreqEst'};
        Control(2).Parameters = {'Process',evParms.proc.procIndCFIIm,...
                                 'srcData','signedColor','persistMethod',...
                                 'dynamic','persistLevel',P(PIndCFI).persfreq};
        Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
        Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
        set(findobj('tag', [slider.perspwr.pos 'Edit']),'String',num2str(P(PIndCFI).persfreq,'%3.0f'));
        set(findobj('tag', slider.perspwr.tag),'Value',P(PIndCFI).persfreq);
        evParms.state.DopState = 'freq';
        % Set modified Process attributes in base Matlab environment.
        Process(evParms.proc.procIndCFI).method = 'computeCFIFreqEst';
        for k = 1:2:length(Process(evParms.proc.procIndCFIIm).Parameters)
            if strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'srcData'), ...
                  Process(evParms.proc.procIndCFIIm).Parameters{k+1} = 'signedColor';
            elseif strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'persistMethod'), ...
                  Process(evParms.proc.procIndCFIIm).Parameters{k+1} = 'dynamic';
            elseif strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'persistLevel'), ...
                  Process(evParms.proc.procIndCFIIm).Parameters{k+1} = P(PIndCFI).persfreq;
            end
        end
    case 2  % Power mode
        newMap = grayscaleCPAmap;
        newMap(1:128,:) = currentMap(1:128,:);
        for k = 1:2:length(Process(evParms.proc.procIndCFIIm).Parameters)
            if strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'persistLevel'), ...
                  P(PIndCFI).persfreq = Process(evParms.proc.procIndCFIIm).Parameters{k+1}; end
        end
        Control(1).Parameters = {'Process',evParms.proc.procIndCFI,'method','computeCFIPowerEst'};
        Control(2).Parameters = {'Process',evParms.proc.procIndCFIIm,...
                                 'srcData','unsignedColor','persistMethod','simple','persistLevel',...
                                 P(PIndCFI).perspwr};
        Control(3).Parameters = {'DisplayWindow',1,'colormap',newMap};
        Control(4).Parameters = {'ImageBuffer',1,'lastFrame',0};
        set(findobj('tag', [slider.perspwr.pos 'Edit']),'String', num2str(P(PIndCFI).perspwr,'%3.0f'));
        set(findobj('tag', slider.perspwr.tag), 'Value', ...
                          P(PIndCFI).perspwr);
        evParms.state.DopState = 'power';
        Process(evParms.proc.procIndCFI).method = 'computeCFIPowerEst';
        for k = 1:2:length(Process(evParms.proc.procIndCFIIm).Parameters)
            if strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'srcData'), ...
                  Process(evParms.proc.procIndCFIIm).Parameters{k+1} = 'unsignedColor';
            elseif strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'persistMethod'), ...
                  Process(evParms.proc.procIndCFIIm).Parameters{k+1} = 'simple';
            elseif strcmp(Process(evParms.proc.procIndCFIIm).Parameters{k},'persistLevel'), ...
                  Process(evParms.proc.procIndCFIIm).Parameters{k+1} = P(PIndCFI).perspwr;
            end
        end
end
evParms.P = P;
assignin('base','evParms', evParms);
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

%CFIDopPowerThreshold - Doppler Power change
uscmsg(execState);
return;
doCFIInternalProc = evalin('base', 'doCFIInternalProc');
if doCFIInternalProc
  Process = evalin('base','Process');
  proc.procIndCFI = evalin('base','proc.procIndCFI');
  for k = 1:2:length(Process(proc.procIndCFI).Parameters)
    if strcmp(Process(proc.procIndCFI).Parameters{k},'pwrThreshold'), Process(proc.procIndCFI).Parameters{k+1} = UIValue; end
  end
  assignin('base','Process',Process);
  % Set Control.Command to set Doppler threshold.
  Control = evalin('base','Control');
  Control.Command = 'set&Run';
  Control.Parameters = {'Process', proc.procIndCFI, 'pwrThreshold',UIValue};
  assignin('base','Control', Control);
else
  assignin('base', 'pwrThresh', UIValue);
  P = evalin('base','evParms.P');
  PIndCFI = evalin('base','PIndCFI');  
  P(PIndCFI).pwrThres = UIValue;
  assignin('base', 'evParms.P', P);  
end
%CFIDopPowerThreshold - Doppler Power change
%qqq
%DopPowerThreshold - Doppler Power change, for SDLarge
evParms = evalin('base', 'evParms');
if evParms.largeSDParms.doSDLargeInternalProc
  Process = evalin('base','Process');
  procIndThis = evParms.proc.procIndSDLarge;
  for k = 1:2:length(Process(procIndThis).Parameters)
    if strcmp(Process(procIndThis).Parameters{k}, 'pwrThreshold')
      Process(procIndThis).Parameters{k+1} = UIValue; 
    end  
  end
  assignin('base','Process',Process);
  % Set Control.Command to set Doppler threshold.
  Control = evalin('base','Control');
  Control.Command = 'set&Run';
  Control.Parameters = {'Process', proc.procIndCFI, 'pwrThreshold',UIValue};
  assignin('base','Control', Control);
  disp(['Set overlay Doppler threshold to: ' num2str(UIValue)]);
else
  evParms.largeSDParms.pwrThres = UIValue;
  assignin('base', 'evParms', UIValue);
  %P = evalin('base','evParms.P');
  %PIndCFI = evalin('base','PIndCFI');  
  %P(PIndCFI).pwrThres = UIValue;
  %assignin('base', 'evParms.P', P);  
end
%DopPowerThreshold

%ColorPriorityLevel - Color Priority change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
proc.procIndCFIIm = evalin('base','proc.procIndCFIIm');
for k = 1:2:length(Process(proc.procIndCFIIm).Parameters)
    if strcmp(Process(proc.procIndCFIIm).Parameters{k},'threshold'), ...
              Process(proc.procIndCFIIm).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.threshold.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process', proc.procIndCFIIm,'threshold',UIValue};
assignin('base','Control', Control);
return
%ColorPriorityLevel

%ColorPersistenceLevel - Color Persistence change
% Set the value in the Process structure for use in cineloop playback.
Process = evalin('base','Process');
proc.procIndCFIIm = evalin('base','proc.procIndCFIIm');
for k = 1:2:length(Process(proc.procIndCFIIm).Parameters)
    if strcmp(Process(proc.procIndCFIIm).Parameters{k},'persistLevel'), ...
              Process(proc.procIndCFIIm).Parameters{k+1} = UIValue; end
end
assignin('base','Process',Process);
% Set Control.Command to set Image.persistLevel.
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process',proc.procIndCFIIm,'persistLevel',UIValue};
assignin('base','Control', Control);

% If PTool window is open, adjust persistLevel1 in Process(proc.procIndCFIIm)
hPTool = findobj('tag','ProcessTool');
if ishandle(hPTool),
    hPNum = findobj('tag','processNum');
    if isequal(get(findobj('tag','processNum'),'Value'),proc.procIndCFIIm)
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

%PRFChange - CFI PRF
evParms = evalin('base','evParms');
P = evParms.P;
PIndCFI = evParms.rcvOut.PIndCFI;
proc.procIndCFI = evalin('base','evParms.proc.procIndCFI');
seqControlIndSetShortTTNACFI = evParms.seq.seqControlIndSetShortTTNACFI;
P(PIndCFI).dopPRF = UIValue;
Trans = evalin('base','Trans');
Process = evalin('base','Process');
demodFreq = evParms.ev.demodFreqCFI;
SeqControl = evalin('base', 'SeqControl');
m = P(PIndCFI).dopNumRays;
SeqControl(seqControlIndSetShortTTNACFI).argument = round(1/(P(PIndCFI).dopPRF*1e-06));

%---------------check Doppler PRF--------------------
currentDepth = evalin('base','Receive(evParms.rcvOut.rcvCFIStartInd(1)).endDepth');
tof = ceil(2*currentDepth/demodFreq);
if SeqControl(seqControlIndSetShortTTNACFI).argument < tof
    SeqControl(seqControlIndSetShortTTNACFI).argument = tof;
    P(PIndCFI).dopPRF = round(1/(tof*1e-06));
    SeqControl(seqControlIndSetShortTTNACFI).argument = round(1/(P(PIndCFI).dopPRF*1e-06));
    fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
    fprintf(['"dopPRF" is adjusted to ' num2str(P(PIndCFI).dopPRF) '\n']);
    UI = evalin('base','UI');
    set(UI(evParms.ind.UIIndCFIPRF).handle(2),'Value', P(PIndCFI).dopPRF);
    set(UI(evParms.ind.UIIndCFIPRF).handle(3),'String',num2str(P(PIndCFI).dopPRF));
end
%---------------check Doppler PRF--------------------

for k = 1:2:length(Process(proc.procIndCFI).Parameters)
    if strcmp(Process(proc.procIndCFI).Parameters{k},'prf'), ...
              Process(proc.procIndCFI).Parameters{k+1} = P(PIndCFI).dopPRF; end
end

evParms.P = P;
assignin('base','P',P);
assignin('base','evParms',evParms);
assignin('base','Process',Process);
assignin('base','SeqControl',SeqControl);

Control(1).Command = 'set&Run';
Control(1).Parameters = {'Process',proc.procIndCFI,'prf',P(PIndCFI).dopPRF};
Control(2).Command = 'update&Run';
Control(2).Parameters = {'SeqControl'};
assignin('base','Control',Control);
return
%PRFChange

%SteerAngle - Doppler Angle will change TX beam and PData
uscmsg(execState);
return;
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


%% **** Callback routines used by External function definition (EF) ****

function setCallback(VsType)
    switch VsType
        case 'Verasonics'
            
            vantageWindow = evalin('base','vantageWindow');
            imageViewer = evalin('base','imageViewer');
            DopRegionLim = evalin('base','DopRegionLim');
            mct  = javaObjectEDT( 'com.verasonics.viewer.tools.mouseclicktool.MouseClickTool', vantageWindow, imageViewer(1));
            mcth = javaObjectEDT( handle(mct, 'callbackproperties' ) );
            mct.setToolEnabled(true);
            mct.lowerLimitX = DopRegionLim.x(1);
            mct.upperLimitX = DopRegionLim.x(2);
            mct.lowerLimitY = DopRegionLim.y(1);
            mct.upperLimitY = DopRegionLim.y(2);
            mcth.ClickedEventCallback = @wbdCallback;
            assignin('base','mct',mct);
        case 'Matlab'
            evalin('base','set(Resource.DisplayWindow(1).figureHandle,''WindowButtonDownFcn'',@nicpseqwt_wbdCallback);');
    end
    return
end
