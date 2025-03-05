function varargout = sdfn20240113(varargin)
%% spectralDoppler  (main) Spectral Doppler processor for pulse Doppler data.
% Copyright 2001-2017 Verasonics, Inc.  All world-wide rights and remedies under all intellectual property laws and industrial property laws are reserved.  Verasonics Registered U.S. Patent and Trademark Office.
%
% DESCRIPTION
%     This routine processes pulse data to obtain classical spectral
%     Doppler output.
%       Usage:
%       1. Environment: Asynchronous processing runacq.c paradigm; called
%       with external process/m-file protocol (See Sequence Programming Manual).
%
%       2. >> spectralDoppler('cleanup');
%       Sets this program to uninitialized state and destroys audio driver
%       object.
%
%       3. >> p=spectralDoppler('spectrogramParameters',framePeriod,dopSweepTimeIndex);
%       Returns spectrogram parameters given frame period in seconds and
%       sweep time control number of 1,2,3, or 4 (fastest to slowest).
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
%       Uses portaudio-based audio driver, mexaudiostream.c.
%
% Assumptions, conditions, and behavior:
%   1. This function is invoked at fixed intervals (usually 20 mSec).
%   2. This function draws lines to the spectral display MatLab window at
%   sweep intervals of 10ms, 20ms, 30ms etc. (multiples of invocation).
%   3. Operating Notes:
%   MACI - Leopard, MatLab R2009b(32bit): must run from shell, invoking matlab with -nojvm option
%       Otherwise "drawnow" causes excessive delays and glitching
%   MACI - SnowLeopard, 32bit, Matlab R2009b(32bit):  runs well from MatLab Gui launch
%   MACI - SnowLeopard, 32bit, Matlab R2010a(32bit):
%               drawnow - 20ms/invocation, occasionally 30ms
%               pause(.0001) - 13ms per invocation, occasionally 30 ms
%               no pause or drawnow - no glitches, but not very smooth.
%%    Vista 64 - runs well.
%
% REVISION HISTORY
%       V1.0:  12/18/2009 JAF stripped down from pwrtas_vs.m
%
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -



if existbasefn('vsExit') & evalin('base', 'vsExit') == 1
  return
end

startTime_clock = clock;
tic;
drawnowUpdateFrameTrigger = 5;

verbose = 0;
testDataGen = 0; %set for fake test data gen
% 1: IQ
% 4: debug vars, including timing
captureData = 5; %set to one to capture IQ data columns in
                 %iqsave_xxx.mat file
captureData = 0; %set to one to capture IQ data columns in iqsave_xxx.mat file
%(is reset later in this file)
%.
%audioDataSave = 1;  %save audio data path diagnostics
audioDataSave = 0;  %save audio data path diagnostics
%can define other capture data sets in switch statement

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

persistent HdWhp
persistent k_chunk sgp subFrameInd IQstate
persistent sounddata sounddataRS sounddataLR %make persistent so reallocation time is reduced
persistent rsParams channelizerState
persistent NoccStartThresh
persistent FsDAC dtDAC
persistent wfState rsState lpfState stateSpecPersist statePhasor
persistent bLPF aLPF startedAudioStream
persistent PWdopPRF SpectralScaling
persistent baselineSave dopWFCutoffNormSave dopSweepTimeIndexSave %use to auto generate the update flag
persistent noiseFloorHistory noiseFloor
persistent spectralDoppler_state launchTime captureVarsString
persistent specWindow SdatAll
persistent computerType invokeCount
persistent sounddataRSCell sounddataCell invocationCount    %debug
persistent prevStartTime startTimeHistory
persistent instanceNumber
persistent SDop
persistent execTimeHistory
%%%%%%%%%%%%%%%%%%%%%%%%%%
global semaphoreKey; semaphoreKey = 24537;
global recordDataFileName; 
persistent vtData

if isempty(vtData)
  vtData = getenv('BNS_DATA');
  if isempty(vtData)
    vtData = 'C:/temp';
    if ~exist([vtData '/tmp'], 'dir')
      mkdir([vtData '/tmp']);
    end
  end
end

if isempty(recordDataFileName)
  recordDataPath = fullfile(vtData, 'nicp', 'tmp');
  if ~exist(recordDataPath, 'dir')
    error(['recordDataPath: ' recordDataPath ' does not exist!'])
  end  
  recordDataFileName = fullfile(recordDataPath, 'recordDataFile.mat');
end

global recordData;
global nonWindowsHost;
persistent semaphoreOpenedState; %This has been changed from global to persistent
                                 % each SD instance create / (mainly opens)
                                 % semaphore only once

if isempty(instanceNumber)
  mfile = mfilename;
  instanceNumber = str2num(mfile(end));
  if instanceNumber > 2
    instanceNumber = 1;
    disp('overriding instance number of sdfn');
  end
  
end

evParms = evalin('base', 'evParms');
forceInit = evParms.SD.forceInit(instanceNumber);

if forceInit
  disp(['Forcing init of SD instance ' num2str(instanceNumber)]);
  clear mexaudiostream
  clear audiostream
  HdWhp=[];
  k_chunk=[]; sgp=[]; subFrameInd=[]; IQstate=[];
  sounddata=[]; sounddataRS=[]; sounddataLR=[];
  rsParams=[]; channelizerState=[]; NoccStartThresh=[];
  FsDAC=[]; dtDAC =[];
  wfState=[]; rsState=[]; lpfState=[]; stateSpecPersist=[]; 
  statePhasor = []; bLPF=[]; aLPF=[]; startedAudioStream=[]; ...
  PWdopPRF=[]; SpectralScaling=[]; 
  baselineSave=[]; dopWFCutoffNormSave=[]; dopSweepTimeIndexSave=[];
  noiseFloorHistory=[]; noiseFloor=[];
  spectralDoppler_state=[]; launchTime=[]; captureVarsString=[];
  specWindow=[]; SdatAll=[];
  computerType=[]; 
  invocationCount=[];
  sounddataRSCell=[]; sounddataCell=[]; %invocationCount 
  prevStartTime=[]; startTimeHistory=[];
  SDop=[];
  execTimeHistory=[];
  whos
end

%First time
if isempty(semaphoreOpenedState)
  %disp('semaphore was closed')
  [dum, osType] = unix('echo $OSTYPE');

  if dum==0 & (strcmp(osType(1:6), 'darwin') | strcmp(osType(1:5), ...
                                                      'linux'))
    nonWindowsHost=1;
    semaphoreOpenedState = 0;
  else
    nonWindowsHost=0;
    %Create OR OPEN semaphore - we use the mexfile associated with semaphore.c
    semaphore('create',semaphoreKey,1); %Initial semaphore value 1;
    semaphoreOpenedState = 1;
  end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%

FrameCountMax = 3000;

%disp(['In instanceNumber ' num2str(instanceNumber)]); 


if isempty(execTimeHistory)
  execTimeHistory = zeros(1, FrameCountMax);
end

if isempty(invocationCount)
    invocationCount = 0;
    startTimeHistory = zeros(1, FrameCountMax);
end
invocationCount = invocationCount+1;

startTime = tic;
if ~isempty(prevStartTime),
    invocationInterval = toc(prevStartTime); % *** remove ; to
                                             % display interval
else
    invocationInterval = 0;
end
prevStartTime = startTime;

% to show invocation delay:
%invocationInterval

startTimeHistory(mod(invocationCount,FrameCountMax )+1 ) = invocationInterval;

stStr = ['startTimeHistory' num2str(instanceNumber)];
assignin('base', stStr, startTimeHistory);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%signature parsing
comm = '';
if ischar(varargin{1}),
    %DUAL COMMANDS
    %first argument is a command:
    comm = varargin{1};
    if length(varargin)>1
        %second argument is a command:
        if ischar(varargin{2}),
            comm2 = varargin{2};   %override
        end
    end
else
    %SINGLE COMMAND but in second argument
    if length(varargin)>1
        %second argument is a command:
        if ischar(varargin{2}),
            comm = varargin{2};   %override
        end
    end
end

%execute subfunction specified in command string variable
if ~isempty(comm),
    
    %is comm a subfunction in this file?
    [ddd,fff]=fileparts(which(comm));
    issubfunc=isequal(fff,mfilename);
    if ~issubfunc,
        error('unrecognized command')
    end
    
    [varargout{1:nargout}]=feval(comm,varargin{2:end});
    
    return %end of command processing mode
    
end %signature check

% jsm put this because persistent variables may not be cleared properly

if exist('SDop', 'var') & ...
      isfield(SDop, 'HSpectrogram') & ...
      ~isvalid(SDop.HSpectrogram)
   spectralDoppler_state = 'not.initialized';
end

if isempty(spectralDoppler_state) | forceInit
    spectralDoppler_state = 'not.initialized';
    %disp('init');
    evParms.SD.forceInit(instanceNumber)
    evParms.SD.forceInit(instanceNumber) = 0;
    assignin('base', 'evParms', evParms);
end

if isempty(computerType),
    computerType=computer;
    invokeCount = 0;
end
invokeCount=invokeCount+1;

%
% .. . . . . . . . . .

wkSpace = 'caller';
scalePortAudio = 10.0; %nominal volume control.
ScaleSpec = 1.0; %adjust spectral
%narrowbanding = evalin(wkSpace,'narrowbanding'); % jsm
NumNoiseHist = 30; %noise floor estimation median window size
narrowbanding = evParms.gate.narrowbanding;

switch nargin,
    
    case 1 % - - - - - - - - -
        %Called by extern object in the script
%        size(varargin{1}) % jsm, IQ buffer size not affected by specific gate
        %get the IQ data:
        IQdataOrigPre = squeeze(varargin{1});
        nPulses = evParms.ev.nDopPRIs; % get no. of PRIs in a
                                       % 'frame'

        % if imag part is zero, matlab will change type to real,
        % which will cause saveiq to try to access a non-existent
        % pointer
        IQdataOrig = IQdataOrigPre(:,:,1:nPulses);        
        
        if isreal(IQdataOrig)
          %disp(['IQ real in instance: ' num2str(instanceNumber)]);
          IQdataOrig = IQdataOrig + j*eps;
        end

%        return % TO DO NOTHING LIVE
        
        IQdata = IQdataOrig;
        nPulses = evParms.ev.nDopPRIs; % get no. of PRIs in a
                                       % 'frame'          
        sizeIQOrig = size(IQdataOrig);

        %need to get because PRF change makes different sizes:       

        SDop.nPulses = nPulses;
        nRange = sizeIQOrig(1);
        
        gate = evParms.gate;
        
        if narrowbanding,
          rangeSubSampleIndex = 1:nRange;
        else         
          rangeSubSampleIndex = gate.rangeSubSampleIndex{instanceNumber};
        end

        if length(sizeIQOrig)==3
          lateralSubSampleIndex = gate.lateralSubSampleIndex{instanceNumber};

          %disp(['nPulses = ' num2str(nPulses)]);
          %disp(['size(IQdata) = ' num2str(size(IQdata))]);

          IQdataSel = IQdata(rangeSubSampleIndex, lateralSubSampleIndex, ...
                             1:nPulses);
          IQdataSel = squeeze(sum(IQdataSel,2));
          NumDepthSummation = length(rangeSubSampleIndex)* ...
              length(lateralSubSampleIndex);            
          
        else
          IQdataSel = IQdata(rangeSubSampleIndex, 1:nPulses);              
          NumDepthSummation = length(rangeSubSampleIndex);
        end

        % can force non-normalization by number of pixels
        %NumDepthSummation = 1;
        
        if length(rangeSubSampleIndex)>1
          IQdata = sum(IQdataSel,1);
        else
          IQdata = IQdataSel;
        end

        % resized now so get new size:
        sizeIQ = size(IQdata);

        if testDataGen ,
            IQdata = generateTestData('tones.clutter.slow',nPulses);
            if rand(1)< 0*.01,
                disp([mfilename,' %%%%%%%%45 '])
                keyboard
            end
        end
        

        if strcmp(spectralDoppler_state,'initialized') % & cond3
            %SDops = evalin(wkSpace, 'SDops');
            %SDop = SDops{instanceNumber};
        else
            SDop = initSDop;
        end

        SDop.noisePersist = evParms.SD.noisePersist;
        SDop.TFrame = evParms.ev.framePeriod_s;
        SDop.PWdopPRF = evParms.ev.dopPRF;
        dopWFCutoffNorm = evParms.SD.dopWFCutoffNorm;
        sweepSpeedIndex = evParms.SD.dopSweepTimeIndex;
        baseLineShiftNorm = evParms.SD.baseLineShiftNorm;
        SDDynRangeDB = evParms.SD.SDDynRangeDB;
        SDop.nIqfft = evParms.SD.NWindow;
        SDop.despeckle = evParms.SD.despeckle;
        
        %set initialization mode flags:
        PWPRFnew = ~isequal(PWdopPRF,SDop.PWdopPRF);
        if PWPRFnew
          disp(['New PRF = ' num2str(SDop.PWdopPRF) ' pulses = ' ...
               num2str(nPulses)]);
        end
        
        PWdopPRF = SDop.PWdopPRF;
        WFnew = ~isequal(dopWFCutoffNorm, dopWFCutoffNormSave);
        dopWFCutoffNormSave =  dopWFCutoffNorm;
        if PWPRFnew
          disp(['New PRF = ' num2str(SDop.PWdopPRF) ' pulses = ' ...
               num2str(nPulses)]);
          %disp(['New dopWFCutoffNorm = ' ...
          %      num2str(dopWFCutoffNorm)]);
%          SDop
        end
        %.
        Sweepnew = ~isequal(sweepSpeedIndex,dopSweepTimeIndexSave);
        dopSweepTimeIndexSave =  sweepSpeedIndex;
        %.
        BLnew = ~isequal(baseLineShiftNorm,baselineSave);
        baselineSave = baseLineShiftNorm;
        
        nuIQRot = baseLineShiftNorm + 0.5;

    otherwise
        error('bad signature switch')
end %nargin switch

instanceNoWithAudioEnabled = evParms.flag.instanceNoWithAudioEnabled;
enableAudio=(instanceNoWithAudioEnabled==instanceNumber);
audioOn = evParms.flag.audioOn & enableAudio;

if ~strcmp(spectralDoppler_state,'initialized')
    
    disp([mfilename,': initializing state.'])

    disp(['In instanceNumber ' num2str(instanceNumber) ...
      ' enableAudio = ' num2str(enableAudio)]);

    PWPRFnew = 1; %override
    
    [sounddataRSCell, sounddataCell] = deal({}); %DEBUG
    
    if enableAudio
      clear mexaudiostream
      clear audiostream
    end
    
    spectralDoppler_state = 'initialized';
    launchTime = startTime_clock;
end

%disable reads/writes from DAC fifo (for testing in VDAS simulation mode).
disableDACDriver = enableAudio~=1 ; %set for no audio
disableSpectralDisplay  = SDop.displayOn~=1; %set for no display

FrameRate = 1.0./SDop.TFrame;
displayTime = NaN;
Nwind = SDop.nIqfft;
  
%% --- Initialize the sweeping display window and Wall filter
% First time through, or any time PRF changes
if PWPRFnew | Sweepnew | BLnew | WFnew,  

  if verbose,
    fprintf('\n')
  end
    
  k_chunk = 0;
  if PWPRFnew ,
    rsState = [];
  end
    
  %calc spectrogram parameters:
  sgp=spectrogramParameters(SDop.TFrame,sweepSpeedIndex);
  THist = sgp.THist;
    
  sgp.nDopPulses = evParms.ev.nDopPRIs; % get no. of PRIs in a 'frame'
                                        % partition the data chunk
                                        % into subrames, one for each line:

  if isnan(sgp.R)
    keyboard
  end
  
  subFrameInd = subframeInd100(sgp,PWdopPRF); %get subframe
                                              %indices within the
                                              %iq chunk

  stateSpecPersist =[];
  if Sweepnew,
    specWindow=[];
  end
    
  [noiseFloorHistory,noiseFloor]=deal([]);
    
  if verbose,
    fprintf('I')
  end

  if enableAudio,
    FsDAC = mexaudiostream(908); % get output sample rate for audio
    dtDAC = 1/FsDAC;
    
    DACFrameSize = mexaudiostream(920); %size of dac frame in stereo samples
    TDACFrame = DACFrameSize * dtDAC ;
    DACFramesPerInvoke =  SDop.TFrame/TDACFrame ; %non integer value
                                                  %two floats per audio sample:
    % how full we want the fifo.
    NoccStartThresh = 1.0*ceil(DACFramesPerInvoke)*DACFrameSize; 
    audioIsOpen = mexaudiostream(3);
    
    disp(['In instanceNumber ' num2str(instanceNumber) ...
      ' audioIsOpen = ' num2str(audioIsOpen)]);
    
    if ~audioIsOpen,
      disp('Audio not yet open, Will now open audio stream:')
      audiostream('open')
     disp('audiostream open call returned.')
    end
  end
  startedAudioStream = 0;
  
  [  channelizerState , wfState , rsState ] = deal([]);
  
  NumDataFrames =  sgp.NumDataFrames;
  nLines =  sgp.NumLines; %actual interpolated lines drawn
    
  %init spectrogram history if size has changed
  %*** experimental force clear
  if ~isequal(size(SdatAll),[SDop.nfft , nLines]),
    SdatAll = zeros( SDop.nfft , nLines);
  end
    
  % make a new figure for this instance
  figureName = ['sdop_figure' num2str(instanceNumber)];
  hFig = findobj('Tag', figureName);
        
  if ~isempty(hFig), close(hFig), end
  TXFreq_MHz = evalin('base', 'Trans.frequency');
  speedOfSound_mps = evalin('base', 'Resource.Parameters.speedOfSound');
  [hFig,hAxs,HImg] = sdop_gui_lite(THist, nLines, SDop.nfft, PWdopPRF, ...
                                   baseLineShiftNorm ,TXFreq_MHz,...
                                   speedOfSound_mps, instanceNumber );
  SDop.hFigure = hFig;
  SDop.hAxes = hAxs;
  SDop.HSpectrogram = HImg;
  if 1
  figPos = figpositionfn;
  if ~isempty(figPos)
    set(hFig,'NumberTitle', 'off');  
    switch evParms.state.CDIState
      case {1,3}
        titStr = evParms.SD.titleFlow;
        rectColorSD1 = evParms.SD.markerColorFlow(1,:);
        rectColorSD2 = evParms.SD.markerColorFlow(2,:);
      case {2,4}
        titStr = evParms.SD.titlePower;      
        rectColorSD1 = evParms.SD.markerColorPower(1,:);
        rectColorSD2 = evParms.SD.markerColorPower(2,:);
    end
    
    setwindowbnsiconfn(hFig);
    evParms = evalin('base', 'evParms', evParms);
    switch instanceNumber
     case 1
      set(hFig, 'position',  figPos.sd2);
      set(hFig, 'name', titStr{1});
      evParms.UI.handleSDFig1 = hFig;
      rh=annotation('rectangle',[.9425 .5 .02 .05], ...
                    'color', rectColorSD1, ...
                    'linewidth', 15);
      evParms.UI.handleSDFig1_rect = rh;
     case 2
      set(hFig, 'position',  figPos.sd1);
      set(hFig, 'name', titStr{2});
      set(gca, 'xtick', [])
      xlabel(gca, '');
      evParms.UI.handleSDFig2 = hFig;
      rh=annotation('rectangle',[.9425 .5 .02 .05], ...
                    'color', rectColorSD2, ...
                    'linewidth', 15);  
      evParms.UI.handleSDFig2_rect = rh;
    end
    assignin('base', 'evParms', evParms);
  end
  end
  
  %sliding FFT window buffer: make sure it's big enough
  % corresponds to 8000 max prf
  [IQstate]=circbuff('initialize', SDop.nfft + 160 ); 
    
  SDop.isInit = 1;
    
  saveIQWithTimeStamp(zeros(size(IQdataOrig)), ...
    evParms.ev.numSaveFramesPerBatch, instanceNumber, ...
                      semaphoreOpenedState);
  %init
  
  switch captureData, %initialize:
   case 0
    %do nothing
   case 1
    IQsave(zeros(nPulses,1),500);%init
   case 2
    captureVars = {...
        'invokeTime','runningTime','displayTime', ...
        'PWdopPRF' ...
                  };
    captureVarsString = vars2VecEvalStr(captureVars);
    IQsave(zeros(length(captureVars),1),100,captureVars);%init
    
   case 3
    captureVars = {...
        'invokeTime','runningTime','displayTime', ...
        'PWdopPRF',  ...
                  };
    captureVarsString = vars2VecEvalStr(captureVars);
    IQsave(zeros(length(captureVars),1),256,captureVars); %init
   case 4
    % jsm added invocation interval
    captureVars = {...
        'invokeTime','runningTime','displayTime', ...
        'drawnowTime' , 'PWdopPRF', 'invocationInterval' ...
                  };
    captureVarsString = vars2VecEvalStr(captureVars);
    IQsave(zeros(length(captureVars),1),256,captureVars);%init
   case 5 % jsm: both IQ and captureVars
    captureVars = {'IQdata(:)', ...
                   'invokeTime','runningTime','displayTime', ...
                   'drawnowTime' , 'PWdopPRF', 'invocationInterval' ...
                  };
    captureVarsString = vars2VecEvalStr(captureVars);            
    ZMat = [zeros(nPulses,1); ...
            zeros(length(captureVars)-1,1)];          
    IQsave(ZMat, 500, captureVars);%init
   otherwise
    disp([mfilename,':error: bad capture option'])
  end
  
  audMute = muter(FrameRate,1); %mute audio
  
else %initialize
    audMute = muter(FrameRate);
end %initialize
IQSlideBlock = zeros(Nwind,sgp.R);
% --- High pass WALL filter design

if WFnew,
    [fproto]= getFilterPrototype('wallfilter');
    [HdWhp.Numerator,HdWhp.Denominator]=filterTranslate(...
        fproto.tf.b,fproto.tf.a,fproto.nuNominal, dopWFCutoffNorm,'high');
end
if ~isequal(size(specWindow),[Nwind,sgp.R]),
    w = hammingvs(Nwind);
    w = w/sum(w);
    specWindow = w(:)*ones(1,sgp.R);
end

SDProcGainDB = 9*log10(Nwind*NumDepthSummation);
baseLineBin = min(SDop.nfft,max(1,round(nuIQRot*SDop.nfft)+1));

%higher is more compression (low level gain).
audCompressionFactor= SDop.audCompressionFactor ;
%spectral display compression method
SDCompressionMethod=SDop.SDCompressionMethod ;

if audioOn,    
    %Update resampler params to PRF= PRFNominal*adjustment
    %until fifo accupancy is ok.
    %Then update resampler params to PRF = PRFEst.
    %%%% fifoState = audiostream('fifostate');
    
    %initialize resampler state:
    if isempty(rsState),
        rsParams = resamplingParams(PWdopPRF,FsDAC);
        if verbose,
            disp([mfilename,':initializing resampler state.'])
        end
        %reasonable (but not optimal) analog resampler design.
        analogResamplerDesign = 'rc.23'; 
        rsState = rs_irrat('initialize', ...
            PWdopPRF , rsParams.FsIntermediate ,nPulses,analogResamplerDesign);
        %second-stage resampler (integer upsample) reconstruction filter
        [fproto]= getFilterPrototype('audio_resampler');
        [bLPF,aLPF]=filterTranslate(...
            fproto.tf.b,fproto.tf.a,fproto.nuNominal,...
            1/(2*rsParams.Qupsample),'low');
        %%%% [bLPF,aLPF] = ellip(5,1.5,60,2/(2*rsParams.Qupsample));
        
    end %init resampler state
    
end% if audio on

audScale = SDop.audioGain * audMute * scalePortAudio;

%finished initialization - - - - - -  - - - - - - - - - - - - -

%SECTION 1: WALL FILTERING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- High pass filter the entire frame of data (Wall Filter)

[IQ,wfState] = filter(HdWhp.Numerator, HdWhp.Denominator, IQdata,wfState);

sMax=6;
%IQ = hankelsvdfiltfn(IQdata, sMax);

if 0
ND = length(IQdata);
P = ceil(ND/2);
A = hankel(IQdata(1:P), IQdata(P:ND));

[U,S,V] = svd(A);

SMod = S;

for s = 1:7
  SMod(s,s)=0;
end

ARec = U*SMod*V';

IQRec = [ARec(:,1); ARec(end, 2:end).'];

IQ = IQRec.'; % .* exp(20 * j*2*pi*(0:ND-1)/ND);

end

%keyboard

if 0
R = A*A';

[U,S,V] = svd(R);

SMod = S;
SMod(1,1)=0;
SMod(2,2)=0;

ARec = U*sqrt(S)*U';

cl = U*SMod*V';


cl = U*U'*IQP(:);



IQP = IQdata(1:P);
clutter = U(:,1)*U(:,1)'*IQP(:) +  U(:,2)*U(:,2)'*IQP(:);

IQP2 = IQdata(P+1:end);
clutter2 = U(:,1)*U(:,1)'*IQP2(:) +  U(:,2)*U(:,2)'*IQP2(:);

IQFilt = [IQP(:); IQP2(:)] - [clutter; clutter2];

%IQFilt = IQP(:) - clutter;
%IQ=IQFilt;
end

%keyboard

%SECTION 2: AUDIO %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% a. channelize doppler in to pos and neg freqs.
% b. convert sample rate by arbitrary fractional amount to nearest rate
% related rationally to DAC sample rate.
% c. upsample by rational fraction rate to DAC sample rate.

% send sound samples for one frame to audio driver

if audioOn,
    %%%%%%%%%%%%%%%%%%
    % --- channelizer operations
    [sounddata,channelizerState]= doppler_channelize(IQ,channelizerState);
    
    % --- convert PRF to next higher freq. that divides the DAC
    % sample playback rate
    [sounddataRS ,rsState] = rs_irrat('resample', sounddata.', rsState);
    
    if audioDataSave,
        sounddataRSCell{end+1}=sounddataRS; %debug
        sounddataCell{end+1}=sounddata;
        if invocationCount==FrameCountMax
            disp([mfilename,': debug save RS audio: audiors_debug.mat'])
            fifodata=mexaudiostream(901,15000);
            fmt=30; %datetime format
            timeAtSave = now;
            vtHome = getenv('BNS_HOME');
            vtData = [getenv('BNS_DATA') '/nicp/'];
            matOutPath = vtData;
            dateFileName = [ matOutPath 'audiors_debug_', ...
                             datestr(now,fmt) '.mat'];
            
            save(dateFileName, 'sounddataRSCell','sounddataCell','fifodata','rsState')
            lslrt(dateFileName);
            invocationCount = 0;
            [sounddataRSCell, sounddataCell] = deal({}); %DEBUG
        end
    end
    
    % --- upsample and LPF by the integer amount to give the DAC sample rate:
    sounddataRS = upsamplevs(sounddataRS,rsParams.Qupsample);
    
    [sounddataRS,lpfState] = filter(audScale*bLPF,aLPF, ...
                                    sounddataRS,lpfState);
    
    %audio compression
    sounddataLR = limiter('soft', ...
        [real(sounddataRS(:)),imag(sounddataRS(:))] , ...
        audCompressionFactor );
    
    %send data to DAC interface fifo:
    if disableDACDriver,
        %for testing on non-realtime environment
        disp([mfilename,':DEBUG: disabled DAC fifo write.'])
    else
        audiostream('load' , sounddataLR);
    end
    
    if ~startedAudioStream,
        if disableDACDriver,
            %for testing on non-realtime environment
            disp([mfilename,': DAC start call disabled.'])
        else
            %check size to start:
            Nocc = 0.5*mexaudiostream(903); %in stereo samples
            if Nocc>=NoccStartThresh ,
                if verbose,
                    %disp([mfilename,': starting audio stream.'])
                    fprintf('S')
                end
                Ncap = (1/2)*mexaudiostream(905); %capacity in stereo samples
                %fraction of fifo occupancy needed to start audio
                startFifoFraction=NoccStartThresh/Ncap ; 
                audiostream('start', startFifoFraction );
                startedAudioStream = 1;
            end
        end
    end

end %if audioOn

% SECTION 3: SPECTRAL DISPLAY %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~disableSpectralDisplay,
    
    k_chunk = rem(k_chunk,sgp.NumDataFrames)+1;
    
    %rotate the frame of data
    [fazor,statePhasor]= phasor( nPulses, nuIQRot, statePhasor);
    % note that IQ should be row vector
    IQrotated = IQ.*fazor;
 %   nuIQRot
    %keyboard

    %spectral processing
    % jsm: how this works
    % 1. next subframe is extracted from this frame
    % 2. it is appended to circular buffer
    % 3. the last Nwind samples of the buffer are retrieved.
    %    (this will end with the current subframe)
    % In this way, overlapping time series are generated 
    for r=1:sgp.R, %collect enough data to draw lines for input chunk
        %R is the spectral slowtime interpolation factor
        % adjust the baseline, using a circular shift
        IQseg = IQrotated(subFrameInd{r});
        %write a subsegment to the buffer
        %(load buffer structure with newest iq data subchunk):
        IQstate = circbuff('write',IQstate,IQseg(:).');
        %get most recent iq data from buffer
        IQsegSlide = circbuff('read',IQstate,Nwind);
        IQSlideBlock(:,r)=IQsegSlide(:); %load multihistory matrix
                                         %for FFT processing
    end
    %fft on all window slides:
    sfBlock = abs(fft(IQSlideBlock.*specWindow,SDop.nfft)).^2;

    % flip direction of flow
    if isequal(SDop.revFlowDir, 1)
        sfBlock = flipud(sfBlock);
    end

    %normalizing spec display data:
    %calc noise floor
    noiseFloorHistory = [median(sfBlock(:,end),1),noiseFloorHistory];
    noiseFloorHistory(NumNoiseHist:end)=[];
    nfhMid = median(noiseFloorHistory);
    if isempty(noiseFloor),
        noiseFloor = nfhMid;
    else
        noiseFloor = noiseFloor*SDop.noisePersist + ...
            nfhMid*(1-SDop.noisePersist);
    end
    
    sfBlock = sfBlock.'/noiseFloor; %normalize to noise floor

    %spectral persistence parameters:
    aSpecPers = [1, -SDop.specPersist];
    bSpecPers = 1-SDop.specPersist;
    %spectral persistence, each bin over slow time:
    [sfBlock,stateSpecPersist] = filter(bSpecPers,aSpecPers, ...
                                        sfBlock,stateSpecPersist) ;
    
    % compress the spectral power
    switch lower(SDCompressionMethod),
      case {'db'}
        noiseFloorColorIndex = -20;
        sfBlock = 10*log10(sfBlock);
        dynRangeComp = SDDynRangeDB+SDProcGainDB;
      case 'power'
        noiseFloorColorIndex = 1;
        sfBlock = sfBlock .^ SDop.compression;
        dynRangeComp = (10^((SDDynRangeDB+SDProcGainDB)/10)).^SDop.compression;
      otherwise
        error('bad switch ')
    end
    %hole-filling processing:
    sfBlock = holefill(sfBlock,SDop.despeckle^(Nwind/SDop.nfft));
    
    %transform to colormap space:
    Sxx = sfBlock.'/dynRangeComp*(SDop.cMapLen - noiseFloorColorIndex)+...
          noiseFloorColorIndex;
    
    %clip to colormap:
    Sxx = min(SDop.cMapLen,max(1,Sxx));
    
    %write baseline:
    Sxx(baseLineBin,:) = SDop.cMapLen/2;
    
    %Sxx(maxFbgn:maxFend) = SDop.cMapLen; %trace of max Frequency in flow
    if ~disableSpectralDisplay,
        imUpdate= paintImage( ...
            'get.line.indices' , SdatAll, Sxx , k_chunk );
        SdatAll(:,imUpdate.columnIndicesRaw) = Sxx;
        if ~isempty(imUpdate.cursorData),
            SdatAll(:,imUpdate.columnIndicesCursor) = imUpdate.cursorData;
        end
        try
        set(SDop.HSpectrogram,'CData', SdatAll); %paint entire spectrogram
        catch
	  disp('Unable to paint spectrogram data!');
          keyboard
        end
        
        tempT = clock;
        switch computerType
            case 'MACI'
                
                if invokeCount==drawnowUpdateFrameTrigger ,
                    pause(.0001) %update gui
                    invokeCount=0;
                end
            case 'MACI64'
                if invokeCount==drawnowUpdateFrameTrigger ,
                    pause(.0001) %update gui
                    invokeCount=0;
                end
            otherwise
                drawnow %update gui
        end
        drawnowTime = etime(clock,tempT);
        
    end
    
end %  spectral display code

if ~evParms.flag.offline 
  saveIQWithTimeStamp(IQdataOrig,evParms.ev.numSaveFramesPerBatch, ...
                      instanceNumber,semaphoreOpenedState);
end

invokeTime = toc;

execTimeHistory(mod(invocationCount,FrameCountMax )+1 ) = invokeTime;
exStr = ['execTimeHistory' num2str(instanceNumber)];
assignin('base', exStr, execTimeHistory);

runningTime = etime(clock,launchTime);

if nargout>0,
    varargout{1} = invokeTime;
end


%diagnostics:
switch captureData,
    case 1
        iqdatasavename = IQsave(IQdata(:));
    case 2
        iqdatasavename = IQsave(eval(captureVarsString),100);
    case 3
        iqdatasavename = IQsave(eval(captureVarsString),256);
    case 4
        iqdatasavename = IQsave(eval(captureVarsString),256);
    case 5
        iqdatasavename = IQsave(eval(captureVarsString),500);
    otherwise
        iqdatasavename='';
end

if ~isempty(iqdatasavename),
    disp([mfilename,':saved diagnostics in file: ',iqdatasavename ...
         ])
end

%disp('called');
%invokeTime



return

end %main
%% End of spectralDoppler.m

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% functions % functions % functions % functions %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a = muter(FrameRate,reset);

persistent a0
if isequal(FrameRate,'cleanup')
    [a0]=deal([]);
    return
end

%muter: mute audio
if nargin<2,
    reset =[];
end
if isempty(reset),
    reset = 0;
end

if isempty(a0),
    a0 = 1.0;
end

alpha = exp(-5/FrameRate);

if reset,
    a0 = 1.0;
    a = 0.0;
    return
end

a0 = a0*alpha;

a = 1.0 - a0;

end %muter


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function s=vars2VecEvalStr(varsList);
%used by data capture utility
s=strcat(varsList,';');
s=[ '[',cat(2,s{:}),']' ];
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Xf=holefill(X,alpha);
%spectral hole filler - set to zero for no holefill
% X must be a stack of row vectors
Xf = zeros(size(X));
Nr = size(X,1);
%note: envdet requires "forgetting factor"
for k=1:Nr,
    [y1,y2]=envdet(1-alpha,X(k,:));
    Xf(k,:) = min([y1;y2]);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function iiCell = subframeInd100(sgp,prf);
%subframeInd100
%partition the data chunk into subframes, one for each line:
%splits up data indices in a frame interval
%R = interpolation factor
R=sgp.R;
TFrame = sgp.TFrame;

% jsm: prf*TFrame is the full number of PRIs in frame
% if CFI is running, this needs to be reduced

%numPulses = prf*TFrame; % original way, rephrased jsm
numPulses = sgp.nDopPulses; % this is actual number of pulses
%nPulses = evalin('base','nDopPRIs'); % get no. of PRIs in a 'frame'

ii = round(linspace(1/R,1,R)'*numPulses); % end indices
II = [[1;ii(1:end-1)+1],ii]; % start indices
for k=1:R,
  try
    iiCell{k} = II(k,1):II(k,2);
  catch 
    keyboard
  end
  
    
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function varargout=paintImage(method,varargin);

persistent Zn2

if isequal(method,'cleanup')
    [Zn2]=deal([]);
    return
end

switch lower(method)
    
    case {'get.line.indices'}%method %%%%%%%%%%%%%%%%%
        [ SdatAll,Sdat,lineIndex]=deal(varargin{:});
        assignin('base', 'SdatAll', SdatAll);        
        Nfft = size(Sdat,1);
        
        if size(Zn2,1)~=Nfft,
            Zn2 =     ones(Nfft,1)*[255 0];
        end
        %convert frame to line0
        R = size(Sdat,2); %num lines per data frame
        NumDataFrames = size(SdatAll,2)/R;
        klineBase = (lineIndex-1)*R;
        klineIndStart=klineBase+ [1:R];
        
        %draw cursor lines:
        cursorIndices = [];
        cursorData = [];
        if lineIndex<NumDataFrames-2,
            cursorIndices = (klineIndStart(end)+[1:2]);
            cursorData = repmat(Sdat(:,R)*2,1,2);
        end
        
        %collect image update info:
        imUpdate.columnIndicesRaw = klineIndStart ;
        imUpdate.columnIndicesCursor = cursorIndices ;
        imUpdate.cursorData = cursorData;
        
        varargout{1} = imUpdate;
        
    case {'set.image.data.local'}%method %%%%%%%%%%%%%%%%%%%%%%%%%%
        %note: only for illustration.  Do this by inline code:
        entireSpectrogramData = varargin{1};
        newSpecData = varargin{2};
        imUpdate = varargin{3};
        
        SdatAll = entireSpectrogramData ;
        SdatAll(:,imUpdate.columnIndicesRaw) = newSpecData;
        SdatAll(:,imUpdate.columnIndicesCursor) = imUpdate.cursorData;
        
        varargout{1} = SdatAll;
        
    case {'set.image'}%method %%%%%%%%%%%%%%%%%%%%%%%%%%
        
        imageHandle = varargin{2};
        
        set(imageHandle,'CData',varargin{1});
        
    otherwise %method %%%%%%%%%%%%%%%%%%%%%
        error('bad switch')
end %method switch %%%%%%%%%%%%%%%%%%%%%%%%
end %paintimage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hFig, hAxs, HImg] = sdop_gui_lite(varargin)
%sdop_gui_lite Spectral Doppler GUI
%.
% [hFig, hAxs, HImg] = sdop_gui_lite(tSweep, nLine, nFreq,
% pwPrf,<baselineshift>, instanceNumber)

m2cm = 100;MHz2Hz = 1e6;


yUnits = 'frequency';
[tSweep, nLine, nFreq, pwPrf, baselineshift, TXFreq_MHz,speedOfSound_mps,instanceNumber]=deal(varargin{:});

% Specified parameters
pX = (600/nLine)*1.25; %pixels per x-unit
pY = 1.0; %pixels per y-unit
nX = nLine;
nY = nFreq;

% Y-axis is scaled from -nY/2+1:nY/2 and labels are in units of velocity
freq2vel = m2cm*speedOfSound_mps/(2*TXFreq_MHz*MHz2Hz); %frequency (Hz) to velocity (cm/s) scale factor
YLim = [-1/2 1/2]*nY;
nYTick = 9; % must be odd

% Colormapping
cMapLen = 128;

CMap = bone(cMapLen );

hAxs = [];

%
if isempty(baselineshift),
    baselineshift = 0.0;
end

hFig = figure;
figureName = ['sdop_figure' num2str(instanceNumber)];
set(hFig,'Toolbar','none','MenuBar','none', ...
    'Position',[100 20 pX*nX+100 pY*nY+100],...
    'Tag',figureName );
switch yUnits,
    case 'frequency'
        yvec = (([0:nY-1]-(nY/2))/nY - baselineshift)* pwPrf;
        yLabelStr = 'Doppler Frequency (Hz)';
    case 'velocity'
        yvec = (([0:nY-1]-(nY/2))/nY - baselineshift)* pwPrf*freq2vel;
        yLabelStr = 'Doppler Velocity (cm/s)';
    otherwise
        error('bad switch')
end
tvec = linspace(0,tSweep,nX);
imageHandle = image(tvec,yvec,zeros(nY,nX));

ylabel(yLabelStr)
xlabel('time (s)')
set(gca,'YDir','normal');
colormap(CMap);

HImg = imageHandle;

drawnow('expose');

end%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sgp=spectrogramParameters(TFrame,sweepTimeInd);
% spectrogramParameters: helper fcn for pwrtas_vs and spectral doppler script
% sgp=spectrogramParameters(TFrame,sweepTimeInd);
% Example:
% sgp=spectrogramParameters(.02,1);
% sgp
%
% sgp =
%
%            TFrame: 0.0200
%        THistAllow: [2 3 4 6]
%          NumLines: 600
%             THist: 2
%     NumDataFrames: 100
%                 R: 6

%jaf 11/1/2009

%round to nearest microsecond:

sgp.TFrame = TFrame;

TFrame_uSec =  round(TFrame*1e6);

%sweeptime parameters are adapted here:
% Want parameters to satisfy the following restriction:
% NumLines/NumDataFrames = NumLines/(THist/TFrame) = integer only

switch TFrame_uSec ,
    
    case 10000 ,
        sgp.THistAllow = [2 3  6];
        sgp.NumLines = 600;
        
    case { 20000 },
        sgp.THistAllow = [2 3 4 6];
        sgp.NumLines = 600;
        
    case { 30000 },
        sgp.THistAllow = [2  3   6];
        sgp.NumLines = 600;
        
    case { 40000 },
        sgp.THistAllow = [ 2 3 4 6 ];
        sgp.NumLines = 600;
        
    case { 50000 },
        sgp.THistAllow = [ 2 3 5 6 ];
        sgp.NumLines = 600;
        
    otherwise
        THAnom_uSec=[2:6]*1e6;
        ndfnom=THAnom_uSec./TFrame_uSec;
        ndf = round(ndfnom);
        LPSnom = 600;
        R =  round(LPSnom./ndfnom) ;
        LPS = R.*ndf;
        THA_uSec = TFrame_uSec.*ndf;
        THA = THA_uSec*1e-6;
        
        sgp.THistAllow = round(THA);
        sgp.NumLines = LPS(sweepTimeInd);


end


THistAllow = sgp.THistAllow;
numAllowedTHist = length(THistAllow);
if  sweepTimeInd > numAllowedTHist,
    error('sweep speed out of range.'),
end
THist = THistAllow(sweepTimeInd);
sgp.THist = THist; %seconds per sgram
NumDataFramesPerSGram = fix(THist/sgp.TFrame); %frames/sgram
sgp.NumDataFrames =  NumDataFramesPerSGram;
sgp.R = round(sgp.NumLines/NumDataFramesPerSGram);  %lines/frame

if ~any(THist==sgp.THistAllow),
    error(['Time history length must be one of: ',num2str(THistAllow)])
end

end %func
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function test;
%test audio and spectrogram with siren IQ data

specDoppFunc = @spectralDoppler; %simplified for delivery

load('L7-4SpectralDoppler') %get parameters

clear mexaudiostream

WavesPerFrame = 6;
NsampPF = framePeriod_s*dopPRF;
iqchunk = getAudioTestDataChunk('siren.1',NsampPF,[],500); %prime
iqchunk=[1 i]*reshape(iqchunk,2,NsampPF);%convert frome stereo to analytic signal
starttime = clock;
schedFrame=0;
for k=1:800, %frame loop
    
    
    %generate test IQ signal:
    iqchunk = getAudioTestDataChunk('siren.1');
    iqchunk=[1 i]*reshape(iqchunk,2,NsampPF);%convert frome stereo to analytic signal
    IQtest= ones(4,1)*iqchunk;
    %add noise and scale:
    IQtest = 8000*(3*IQtest + 1*(randn(size(IQtest)) + i*randn(size(IQtest)))/sqrt(2) );
    
    extime = feval(specDoppFunc,IQtest);
    
    %wait for next sched. frame time before generating next test frame:
    while (k-1)>(schedFrame-15),
        runningTime = etime(clock,starttime);
        schedFrame = runningTime/framePeriod_s;
        switch computer,
            case {'PCWIN','PCWIN64'}
                if rand(1)<.0001,
                    fprintf('.')
                end
            case {'MACI','MACI64'}
                fprintf('.')
            otherwise
                fprintf('.')
        end
    end
    if rem(k,50)==0,disp(['test frame: ',num2str(k)]),end
    
end

clear mexaudiostream


disp('...test done.')
%%%%%%%%%%%%%%%%%%%%

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
function x=getAudioTestDataChunk(testSignalType,varargin);
%getAudioTestDataChunk: returns stereo audio data chunk as row vector
%usage:
% initial call:
% x = getAudioTestDataChunk(testSignalType,Nload,<Nblocks>);
% x ~ [1 x Nload ]
% default Nblocks is  500
% subsequent calls:
% x = getAudioTestDataChunk(testSignalType);
% methods:
% testSignalType = 'siren.1' - continuous
%

% john flynn 11dec2009

persistent PhzVecs kload Nload NFrameCyc vibratoAmount
if isequal(testSignalType,'cleanup')
    [PhzVecs ,kload ,Nload ,NFrameCyc ,vibratoAmount]=deal([]);
    return
end


%parse inputs:
Kvai = length(varargin);kvai = 1; %use kvai and template below:
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; NloadIn = tempVar;
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; NFrameCycIn = tempVar;
if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
else tempVar=[];
end;kvai=kvai+1; vibratoAmountIn = tempVar;


if isempty(PhzVecs)
    disp([mfilename,': initializing audio test signal gen...'])
    
    if isempty(NFrameCycIn),
        NFrameCycIn = 500;
    end
    if isempty(vibratoAmountIn),
        vibratoAmountIn = 2;
    end
    if isempty(NloadIn),
        error('must specify load size')
    else
        Nload = NloadIn; % samples per DAC frame
    end
    NFrameCyc = NFrameCycIn;
    vibratoAmount = vibratoAmountIn;
    
    %siren phase history:
    PhzVecs = 2*pi*3*[0:NFrameCyc*Nload-1]/(Nload)+ ...
        (sin(2*pi*3*[0:NFrameCyc*Nload-1]/(NFrameCyc*Nload))*vibratoAmount)*2*pi;
    PhzVecs = reshape(PhzVecs,Nload,NFrameCyc)';
    
    kload = 0;
    
end

Nblock = size(PhzVecs,1);
switch lower(testSignalType),
    case {'siren.1'}%testSignalType %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kload = mod(kload+1-1,Nblock)+1;
        phzvec=PhzVecs(kload,:);
        x = [cos(phzvec);sin(phzvec)];
        x=x(:).';
    otherwise
        error('bad test signal type name')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [y,stateOut]=doppler_channelize(x,state);
% doppler_channelize: split Doppler IQ into L/R by spectral halves.
% Output is coded into real/imag as left/right
% usage:
%  [y,stateOut]=doppler_channelize(x);
%  [y,stateOut]=doppler_channelize(x,state);
%
% Output    Spectrum
% - - - - - - - - -
% real      positive
% imag      negative
%
% ----------------------------------------------------
% John Flynn 16Feb2009
% (c) 2009 VeraSonics, Inc.

%initialize
ChannelizerBWFactor = 1.0; % range = (0.0 -  1.0)

if ischar(x),
    %special command processing
    command = x; clear x
    switch lower(command)
        case {'test'} %comm switch %%%%%%%%%%%%%%%%%%%%%%%
            
            x = dopplersignal(30000);
            X = reshape(x ,3,10000);
            
            X = conj(X); %switches audio to opposite channel
            
            [Y(1,:),dcs] = doppler_channelize(X(1,:));
            [Y(2,:),dcs] = doppler_channelize(X(2,:),dcs);
            [Y(3,:),dcs] = doppler_channelize(X(3,:),dcs);
            
            y = Y.';
            y = Y(:);
            
            y2 = upsamplevs(y,5);
            
            ylr = [real(y2(:)),imag(y2(:))];
            
            %should have sound in left or right channel only
            soundsc(ylr,44100);
            
        otherwise
            error('bad switch ')
    end %comm switch %%%%%%%%%%%%%%%%%%%%%%%%%
    
    return;
    
end%command check

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% normal processing


if nargin<2,
    state = [];
end

if isempty(state),
    % init state:
    state.fshifters = [] ; %note: all stages are same since all sizes are same
    state.lpfLeft = [];
    state.lpfRight = [];
    
    %    filter design:
    [fproto]= getFilterPrototype('doppler_channelize');
    [b,a]=filterTranslate(...
        fproto.tf.b,fproto.tf.a,fproto.nuNominal,  0.25*ChannelizerBWFactor ,'low');
    
    %%%%[b,a]=ellip(LPFOrd,LPFPassBAtten,LPFStopBAtten,2*(0.25*ChannelizerBWFactor));
    state.b = b;
    state.a = a;
    
end
b = state.b;
a = state.a;

stateOut = state;

[xhi,fshifters] = freqshift(x,-0.25,state.fshifters);
[xlo] = freqshift(x,0.25,state.fshifters);
stateOut.fshifters = fshifters;

[xhif , stateOut.lpfRight]= filter(b,a,xhi,state.lpfRight);
[xlof , stateOut.lpfLeft] = filter(b,a,xlo,state.lpfLeft);
[yr] = freqshift(xhif,0.25,state.fshifters);
[yi] = freqshift(xlof,-0.25,state.fshifters);

y = real(yr)+i*imag(yi);

end%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function xc = limiter(method,x,varargin);
%limiter: saturation appropriate for audio compression/compressor
% compression/limiting is w.r.t output range of +/- 1.0
%usage:
% xc = limiter('soft',x,varargin);
% xc = limiter('hard',x,varargin);

%jaf 14apr2009

if isempty(method)
    method = 'hard';
end

Kvai = length(varargin);kvai = 1; %use kvai and template below:

switch lower(method),
    case {'hard'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xc = x;
        [posLim,negLim]=deal(1,-1);
        
        x(find(x>posLim)) = posLim;
        x(find(x<negLim)) = negLim;
        
    case {'soft'} %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; compressionFactor = tempVar;
        
        if isempty(compressionFactor)
            compressionFactor = 1.0;
        end
        
        xc = atand( x * compressionFactor) /(90 );
        
    otherwise %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        error('bad switch ')
        
end %method %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [y,stateOut,s] = freqshift(x,nu,state);
%freqshift: shift freq of input vector by nu.
%usage:
% [y,stateOut,<shifterVec>] = freqshift(x,nu,state);
% x is a vector of signal samples.
% nu is normalized frequency.
%

%jaf 16feb2009


[m,n]=size(x);
isRow = m==1;

if min(m,n)>1
    error('matrix input not supported')
end
if nargin<3,
    state = [];
end

if isempty(state),
    state.nextStartPhasor = 1;
end

L = length(x);


if isRow,
    s = exp((2*pi*i*nu)*[0:L])*state.nextStartPhasor;
else
    s = exp((2*pi*i*nu)*[0:L]')*state.nextStartPhasor;
end

y = x.*s(1:end-1);

stateOut.nextStartPhasor = s(end);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [x,stateOut]=phasor(N,nu,state);
%phasor: complex phasor with frequency nu, N samples
%usage: [x,state]=phasor(N,nu,state);
%N - number of samples
%nu normalized freq.
%state- init for recurrence

%john flynn 12/7/2009
if nargin<3,
    state = [];
end

if ischar(N),
    comm = N;
    switch lower(comm)
        
        case 'test'
            
            [x1, s]=phasor(12,.134);
            [x2,s]=phasor(6,.134,s);
            [x12]=phasor(18,.134);
            
            dx=x12-[x1,x2];
            err=norm(dx)/norm(x2);
            
            x = err;
            disp([mfilename,':',comm,': test errors: ',num2str(err)])
            
        otherwise
            
            error('bad switch')
            
    end
    
    return
    
end

if isempty(state)
    state = 0.0;
end
%   2*pi*i*nu*[0:N] + 2*pi*i*phi
% = 2*pi*i*(nu*[] + phi)
% = 2*pi*i*(nu*([]+phi/nu))
phz = (2*pi*nu*i)*[(0+state):(N-1+state)] ;
stateOut = N+state;
x=exp(phz);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
function rsParams = resamplingParams(PRF,FsDAC);
% resamplingParams: resampling parameters
%.
% Given PRF and FsDAC, returns a structure with:
% FsIntermediat = Intermediate resample rate;
% Qupsample = upsampling factor
%

%jaf 17feb2009

if FsDAC<PRF,
    error('PRF is higher than specified DAC sample rate')
end

MaxQ = 48;

Qset = 2:MaxQ;
FsDACsub = FsDAC./Qset;
ff = find(FsDACsub>PRF) ;
Index_Fs_leastGT_PRF = ff(end);

rsParams.Qupsample= Qset(Index_Fs_leastGT_PRF);
%upsampler rate from smallest subsample GT PRF

rsParams.FsIntermediate = FsDACsub(Index_Fs_leastGT_PRF);
rsParams.FsDAC = FsDAC;
rsParams.PRF = PRF; %input Fs

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout]=circbuff(method,varargin);
%circbuff: manage circular buffer.
% read does not advance pointers.
%signatures -
% [state]=circbuff('initialize',Nbuff);
% [state]=circbuff('write',state,dataIn);
% [dataOut]=circbuff('read',state,Nread);

%john flynn 11/5/2009
Kvai = length(varargin);kvai = 1; %use kvai and template below:

switch lower(method),
    case {'initialize'} %--------------------------
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; NBuff = tempVar;
        
        state.NBuff = NBuff;
        state.isInitialized = 1;
        state.nextWriteIndex = 1;
        state.data = zeros(1,state.NBuff);
        
        varargout{1}=state;
        
    case {'write'} %--------------------------
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; state = tempVar;
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; dataIn = tempVar;
        
        Nnew = length(dataIn);
        wrInd = [0:Nnew-1]+state.nextWriteIndex;
        wrInd = restrictIndex(state,wrInd);
        
        state.data(wrInd) = dataIn;
        state.nextWriteIndex = restrictIndex(state,state.nextWriteIndex+Nnew);
        state.isInitialized = 0;
        
        varargout{1}=state;
        
    case {'read'} %--------------------------
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; state = tempVar;
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; Nread = tempVar;
        
        if state.isInitialized,error('attempt to read from empty circular buffer '),end;
        
        K=state.nextWriteIndex;
        ri = [K-1-Nread+1:K-1];
        ri = restrictIndex(state,ri);
        
        varargout{1} = state.data(ri);
        
        
    otherwise %--------------------------
        error( 'bad switch' )
end %--------------------------


    function ir=restrictIndex(state,indVec);
        ir = mod(indVec-1,state.NBuff)+1;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dataFileName]=IQsave(x,NframesIn,saveVarNames);
%IQsave: capture IQ data into buffer and save into file.
%.
%initialize:
% IQsave(dataColumnSnapshop,NumHist);
%.
% accumulate and auto save if history is full:
% saveFileName = IQsave(dataColumnSnapshop);
% saveFileName is empty unless save occurred.

%jaf 05apr2009
%

verbose = 1;

persistent X nframes Nframes captureVars

if nargin<3,
    saveVarNames=[];
end

if isequal(x,'cleanup')
    [X ,nframes, Nframes]=deal([]);
    return
end

if isempty(captureVars)
    if isempty(saveVarNames),
        [MX,NX]=size(X);
        saveVarNames = cell(1,MX);
        saveVarNames(:) = {''};
    end
    captureVars = saveVarNames;
end

doInit = length(x)~=size(X,1);
if doInit,
    %initialization
    if nargin<2,
        NframesIn = [];
    end
    if isempty(NframesIn),
        NframesIn = 300;
    end
    Nframes = NframesIn;
    X = single(zeros(length(x),Nframes));
    nframes = 0;
else % jsm put the following in conditional, else
     % X(:,1) will always have only zeros
  nframes = nframes+1;
  nframes = rem(nframes-1,Nframes)+1;
  X(:,nframes) = single(x(:));
end

saveBuffer = nframes==Nframes ;
if saveBuffer,
  disp('Saving buffer');
    fmt=30; %datetime format
    timeAtSave = now;
    vtData = getenv('BNS_DATA');
    matOutPath = [vtData '/nicp/iqdata/'];             
    dataFileName = [matOutPath 'iqsave_',datestr(now,fmt) '.mat'];
    creator = mfilename;
    callingStack = dbstack;
    savelist = {'X' , 'Nframes' , 'creator', 'dataFileName' , 'callingStack' ,'captureVars', 'timeAtSave'};
    save(dataFileName, savelist{:});
    lslrt(dataFileName);

    if verbose,
        disp(' ')
        disp([mfilename,':captured data in file: ', dataFileName,'.mat'])
    end
else
    dataFileName = '';
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function SDop = initSDop;
%defaults
SDop = struct('PWdopPRF', [], ...
    'TFrame' , [] , ...
    'nSampleVolume', 1, ...
    'nPulses', [], ...
    'nfft', 256 , ... %FFT size
    'nIqfft', 35 , ... %FFT time window
    'revFlowDir', 0 , ... %flips axis
    'specPersist', 0.3 , ... %spectral slowtime smoother, 0.0 to 1.0
    'despeckle',0.45 , ...  %hole filler, 0.0 to 1.0
    'noisePersist', 0.95 , ... %noise floor smoother, 0.0 to 1.0
    'displayOn', 1, ...
    'cMapLen', 128, ...
    'audioOn', 1 , ... %set to 0 to disable audio
    'audioGain', 3e-6, ...
    'audCompressionFactor' , 1.5, ... %higher is more compression (low level gain).
    'compression', 0.5, ... %spectral display compression level for power method,
    ... % if set to 1.0 ==> makes the colormap index units be watts into 1 Ohm load
    'SDCompressionMethod' , 'power' , ... %set to 'power' or 'db'
    'isInit', 0);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fproto]= getFilterPrototype(filterPurpose);
%getFilterPrototype: retrieve digital filter prototype for specific
% spectral doppler task.

persistent filterDesigns
if isempty(filterDesigns)
    load spectraldoppler_filtercoeff
end

fproto =filterDesigns.(filterPurpose);

end%
%%%%%%%%%%%%%%%%%%%%%%%%
function filterDesigns = filterDesignSpecs
%this is called by specDoppDesigner, which can be isolated from remainder
%of spectral doppler package for a compilation not dependent on the Signal Processing Toolbox

%call filter design functions with nominal critical freq,
% and desired attenuation params, in lowpass form.
% Will transform as needed to desired critical frequency at
% run-time.
filterDesigns = struct;
nuNominal = 0.25;
style = 'low'; %will be transformed to desired style in runtime.
% - - - -
%custom part:
filterPurpose = 'audio_resampler'; %at top level
callingSignature= { 5,1.5,60, nuNominal*2 ,style};
filterHandle = 'ellip';
%common:
fds.nuNominal = nuNominal;
fds.filterPurpose = filterPurpose;
fds.filterHandle = filterHandle;
fds.callingSignature = callingSignature;
filterDesigns.(filterPurpose) = fds;

% - - - -
%custom part:
filterPurpose = 'doppler_channelize'; %in a subroutine
callingSignature={5,1.5,30, nuNominal*2 ,style};
filterHandle = 'ellip';
%common:
fds.filterPurpose = filterPurpose;
fds.filterHandle = filterHandle;
fds.callingSignature = callingSignature;
filterDesigns.(filterPurpose) = fds;

% - - - -
%custom part:
filterPurpose = 'wallfilter';
callingSignature= { 7, 65  nuNominal*2 ,style };
filterHandle = 'cheby2';
%common:
fds.filterPurpose = filterPurpose;
fds.filterHandle = filterHandle;
fds.callingSignature = callingSignature;
filterDesigns.(filterPurpose) = fds;

varargout{1}= filterDesigns;

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y = generateTestData(method,ppf)
t = [0:(ppf-1)]/ppf  ; %time vector to wrap every frame
iphz = 2*pi*1i*t ;
audioScale = 40000;
clutterLevel = 330;
switch method
    case 'tones.1'
        flow = exp(iphz*20) +  exp(-iphz*18)  ;
        y=audioScale*(  flow );
    case 'tones.clutter.slow'
        flow = exp(iphz*8) +  exp(-iphz*9)  ;
        clutter = clutterLevel*exp(iphz);
        y=audioScale*( clutter +  flow  );
    case 'tones.clutter.dc'
        flow = exp(iphz*20) +  exp(-iphz*18)  ;
        clutter = clutterLevel*1;
        y=audioScale*( clutterLevel*clutter +  flow  );
    otherwise
        error('bad switch')
end
end % generateTestData
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = restrictAcqPRF( varargin );

%parse inputs:
Kvai = length(varargin);kvai = 1; %use kvai and template below:

method = 'prf-constraints.1';

switch method
    case 'legacy'
        %inputs: prfnom,tframenom
        [prfnom,tframenom] =deal(varargin{[2 1]} );
        % parameter in seconds or hz
        %produces even number of acq pulses per frame
        ttnanom = 1/prfnom ;
        ttnanom_microsec = ttnanom*1e6 ;
        ttna_microsec = round(ttnanom_microsec);
        ttna = ttna_microsec/1e6;
        prf = 1/ttna;
        if nargin>1,
            ppffnom = prf*tframenom;
            ppf=round(ppffnom);
            %round to higher even
            ppf = ppf + double(rem(ppf,2)==1);
            tframe = ppf/prf;
        end
        varargout{:} = deal( prf,ttna_microsec,ttnanom_microsec,tframe,ppf );
        
    case 'prf-constraints.1' %-------------------
%signature:  FTNom_uSec, PRFDopNom_Hz
% where:  
%  class(FTNom_uSec) == 'char'  ==> use exactly
%  class(FTNom_uSec) == numeric ==> will quantize

        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; FTNom_uSec = tempVar;
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; PRFDopNom = tempVar;
        if kvai<=Kvai, tempVar=varargin{kvai} ; %%%%%%%%%%%%%%
        else tempVar=[];
        end;kvai=kvai+1; verbose = tempVar;
        
          
        nc=[190:600];  %candidate clocks/PRI
        
        if isempty(FTNom_uSec)
%default
            FTNom_uSec =  '2*2*2*2*3*3*5*7*11' ; % 55440 
        end
        if isempty(verbose)
            verbose =  0 ;
        end
        
        %choose closest allowed Frametime
        %frame time in microseconds
        if ischar(FTNom_uSec),
            %use exact
            FT_uSec = eval(FTNom_uSec);
        else
            %will quantize to "best" nearest
            %best allowing choice of prfs without changing frame time
            FTvec_uSec = [ ...
                eval('2*2*2*2*2*2*2*3*3*5*7') ... 40320
                eval('2*2*2*3*3*5*11*13') ... 51480
                eval('2*2*2*2*3*3*5*7*11') ... 55440
                eval('2*2*2*2*3*3*5*7*12') ... 60480
                ];
            [~,indft]= min(abs(FTNom_uSec-FTvec_uSec));
            FT_uSec = FTvec_uSec(indft);
        end
        
        pphalffnom=(FT_uSec/2)./nc;
        pphalffc=round(pphalffnom);
        iii=find(abs(pphalffnom-pphalffc)./abs(pphalffc)<1e-6);
        pphalff = pphalffc(iii);
        pria=nc(iii)*1e-6;
        
        prfa=1.0./pria;
        prfd = (prfa/2);
        prfdnom=round(prfd);
        
        pulseSched.PRFAcq = prfa;
        pulseSched.PRFDop = prfd;
        pulseSched.maxPRFAcq = max(prfa);
        pulseSched.maxPRFDop = max(prfd);
        pulseSched.minPRFAcq = min(prfa);
        pulseSched.minPRFDop = min(prfd);
        pulseSched.pulsesPerFrameAcq  = round(pphalff*2);
        pulseSched.pulsesPerFrameAcqMax  = max(pulseSched.pulsesPerFrameAcq);
        pulseSched.pulsesPerFrameDop  = round(pphalff);
        pulseSched.pulsesPerFrameDopMax  = max(pulseSched.pulsesPerFrameDop);
        pulseSched.PRIAcq_uSec  = nc(iii);
        pulseSched.PRIDop_uSec  = nc(iii)*2;
        pulseSched.FrameTime_uSec  = round(FT_uSec);

        if ~isempty(PRFDopNom)
            [~,indPRF] = min(abs(pulseSched.PRFDop-PRFDopNom));
        else
            indPRF = [];
        end
        pulseSched.indexPRF = indPRF;
        
        if ~isempty( pulseSched.indexPRF)
             fnames=fieldnames(pulseSched);
             for k=1:length(fnames),
                 varname = fnames{k};
                 var = pulseSched.(varname);
                 L=length(var);
                 if L>1
                     pulseSched.(varname)=var(pulseSched.indexPRF);
                 end
             end
        end
        
        if verbose,
            FT_uSec
            plot(prfdnom,'o'),
            grid on
        end
        varargout{1}= pulseSched;

    otherwise
        error('bad switch')
end %method switch
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function unInitPersistentInCaller
% assign empty to all persistent vars in caller wkspc
w=evalin('caller','whos;');

wp={w.persistent};
wp=cat(2,wp{:});
wn = {w(wp==1).name};

for k=1:length(wn),
    evalin('caller',[wn{k},'= [];']);
end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cleanup
%close audio driver:
audiostream('cleanup')
%need to init persist in caller (not in this function)
evalin('caller', 'unInitPersistentInCaller;')

%uninit. subfunctions with persistent vars:
muter('cleanup')
paintImage('cleanup')
getAudioTestDataChunk('cleanup')
IQsave('cleanup')

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [dataFileName]=saveIQWithTimeStamp(x,NframesIn, instanceNumber,...
                                            semaphoreOpenedState)
  %saveIQWithTimeStamp: capture IQ data  into file as controlled by pump monitor
  % For Init: IQsave(zeros(nPulses,1),60);%init
  % For subsequent iqdatasavename = IQsave(IQdata(:),60);
  evParms = evalin('base', 'evParms');
  verbose = 1;
  dataFileName = '';
  persistent X Timestamp TimestampNum nframes Nframes nframesWatchdog ...
    sessionName sessionNameRec pausedMsgPrint stepNo stepTrialNo stepSegNo;
  
  if isequal(x,'cleanup')
    [X ,nframes, Nframes]=deal([]);
    disp('saveIQWithTimeStamp: cleanup');
    nframesWatchdog=[];
    sessionName = '';
    sessionNameRec = ''; % tracked in saveBuffer conditional
    return
  end

  %Check for initialization
  doInit = ~(size(x,1)==size(X,1) & size(x,2)==size(X,2) & ...
             size(x,3)==size(X,3));
  if doInit,
    disp('saveIQWithTimeStamp: doInit');
    pausedMsgPrint=0;
    %initialization
    if nargin<2,
      NframesIn = [];
    end
    if isempty(NframesIn),
      NframesIn = evParms.ev.numSaveFramesPerBatch;
    end
    Nframes = NframesIn;
    X = zeros([size(x) Nframes], 'single');
    Timestamp =  cell(1,Nframes);
    TimestampNum =  zeros(1,Nframes);        
    for i = 1:Nframes 
      Timestamp{1,i} = '';
    end
    nframes = 0;
    nframesWatchdog=0;
    subjectName = '';
    sessionName = '';
    sessionNameRec = '';    
    stepNo = -1;
    stepTrialNo = -1;    
    stepSegNo = 1;
  else% (Not Initialization) 
    global recordData;
    successFlag = getrecordinginfofn(semaphoreOpenedState);

    % (evParms.flag.doWT | evParms.flag.BLineWT) & ...
%    recordData.stepType
%    recordData.UIModePriority
    if recordData.stepType(1) == 'w' & ...
          ~recordData.UIModePriority
      % wall track
      recordData.recordState = 0; % force finish recording
      switchToWT=1
    else
      switchToWT=0;
    end
    
    % push into SDOnly mode without starting record: reduces delay
    % in switching at start of protocol step
    
    if recordData.stepType(1) == 'd' & ~recordData.UIModePriority    
      forceSwitchtoSDOnly = 1;
    else
      forceSwitchtoSDOnly = 0;
    end
    
    
    if evParms.flag.forceRecord | (1 == successFlag)
      % intention is to begin a new file upon resume
      doRecording =  (1 == recordData.recordState || ...
                      3 == recordData.recordState); % start/resume
      % finish recording if we stop recording for any reason
      finishRecording = (forceSwitchtoSDOnly || (0 == recordData.recordState || ...
                          2 == recordData.recordState)& ...
                          1 == evParms.state.recording(instanceNumber));
      evalin('base', ...
             ['evParms.state.recording(' num2str(instanceNumber) ') = ' ...
              num2str(doRecording) ';']);
    
      % Put things in buffer only in started or resumed state
      if evParms.state.SDOnlyInPlace | ...
            ~evParms.flag.forceHighPRFModeForRecordingIQ             
        if evParms.flag.forceRecord | doRecording | finishRecording    
          if ~strcmp(recordData.sessionName, sessionName) % new
                                                          % session
            disp([mfilename ': Detected new session: ' recordData.sessionName]);
            disp([mfilename ': Old session was: ' sessionName]);
            nframes=0;
            X = zeros([size(x) Nframes], 'single');
            sessionName =  recordData.sessionName;
            stepSegNo=1;
            finishRecording = 0; % don't finish any previous recording
          end
                
          nframes = nframes+1;
          nframes = rem(nframes-1,Nframes)+1;
          X(:,:,:,nframes) = single(x);
          dateNum = datenum(now);
          dateStr = datestr(dateNum, 'yyyymmddHHMMSSFFF'); 
          TimestampNum(nframes) = dateNum;
          saveBuffer = (nframes==Nframes) ||  finishRecording;
        
          nframesWatchdog = nframesWatchdog+1;
          nframesWatchdog = rem(nframesWatchdog-1, ...
                                evParms.ev.watchdogInterval_frames)+1;
          watchdogUpdate = (nframesWatchdog == ...
                            evParms.ev.watchdogInterval_frames);
          if watchdogUpdate & instanceNumber==1
            vtData = getenv('BNS_DATA');
            watchdogFilename = fullfile(vtData, evParms.ev.watchdogFilename);
            fidw = fopen(watchdogFilename, 'w');
            if fidw
              fprintf(fidw, dateStr);
              fclose(fidw);
            else
              disp(['Error: Watchdog file: ' watchdogFilename ...
                    ' could not be written.']);
            end                      
          end
          
          if saveBuffer 
            vtData = getenv('BNS_DATA');     
            % check if we need another file to complete recording
            % this step 
            contFlag = (recordData.stepNo==stepNo) & ...
                       (recordData.stepTrialNo == stepTrialNo) & ...
                        strcmp(sessionNameRec, ...
                               recordData.sessionName);
            sessionNameRec = recordData.sessionName;
            
            if contFlag
              disp('Need another file to complete recording this step.');
              stepSegNo=stepSegNo+1;
            else
              disp('First segment.');
              stepSegNo=1;
            end
            
            outPath =  fullfile(vtData, 'nicp', ...
                                recordData.sessionName);
            
            outFilePrePre = fullfile(outPath, ...
                                  [recordData.subjectName ...
                                '_s' num2str(recordData.stepNo) ... 
                                '_t' ...
                                num2str(recordData.stepTrialNo)]);
            
            if instanceNumber == 1 & evParms.largeSDParms.largeSDSave
              % make available for SDLarge recording
              evalin('base', ['currentOutputFilePrefix = ''' outFilePrePre ...
                              ''';']);
            end
            
            outFilePre = [outFilePrePre ...
                          '_sg' num2str(stepSegNo) ...
                          '_i' num2str(instanceNumber)];  

            
            
            % prevent overwrites
            fileFixStr = '';
            cntFix=0;
            cont=1;
            while cont
              iqFileRe = [outFilePre fileFixStr '.re'];
              iqFileIm = [outFilePre  fileFixStr '.im'];                
              if exist(iqFileRe, 'file') | exist(iqFileIm, 'file')
                disp(['*** Warning: ' iqFileRe ' exists!']);
                cntFix=cntFix+1;
                fileFixStr = ['_fix' num2str(cntFix)]; 
              else
                cont=0;
              end            
            end
            
            stepNo = recordData.stepNo;
            stepTrialNo = recordData.stepTrialNo;
            
            PData = evalin('base', 'PData');                    
            PDataThis = ...
                PData(evParms.gate.PDataIndSDStart-1+instanceNumber);
            calledFn = ['saveiqfn_' ...
                        num2str(instanceNumber)];

            if ~exist(outPath, 'dir')
              disp(['*** Output path does not exist. Making it to ' ...
                    'avoid crash. ***']);
              mkdir(outPath);
            end
            
            evParms.this.mode = 'd';
            evParms.this.TWFreq_MHz = evParms.TW(evParms.ind.TWIndSD).Parameters(1);
             
            feval(calledFn, X(:,:,:,1:nframes), evParms, PDataThis, ...
                  iqFileRe, iqFileIm, ...
                  TimestampNum(1:nframes));
            % reset buffer
            X = zeros([size(x) Nframes], 'single');
            nframes=0;

	    if instanceNumber == 1
	      % save ultrasound image window
	      saveas(...
	        evalin('base','Resource.DisplayWindow(1).figureHandle'), ...
		[outFilePre '.fig']);
	    end
	    
          end % saveBuffer
          
          if 0 & saveBufferMat,
            matOutPath = [getenv('BNS_DATA') '/nicp/'];
            if ~exist(matOutPath, 'dir')
              error(['Data path: ' matOutPath ' does not ' ...
                     'exist']);
            end                    
            
            dataFileName = [matOutPath subjectName '_SDI_' ...
                            num2str(instanceNumber) ...
                            '_IQ_',datestr(now,'yyyymmddHHMMSSFFF') ...
                            '.mat'];
            savelist = {'Nframes' , 'Timestamp', 'X'};
            
            save(dataFileName, savelist{:});
            %Reset buffer
            X = single(zeros(length(x),Nframes));
            for i = 1:Nframes 
              Timestamp{1,i} = '';
            end
            lslrt(dataFileName);
            if verbose,
              disp(' ')
              disp([mfilename,':captured data in file: ', dataFileName])
            end %if verbose
          end % if saveBufferMat
          
        end % evParms.flag.forceRecord | doRecording
      else     
        if (doRecording |  forceSwitchtoSDOnly) & ~evParms.state.SDOnlyInPlace      
          evParms = evalin('base', 'evParms');
          evParms = setoverlaysdpriorityfn(evParms);                  
          evParms.state.SDOnlyInPlace=1;
          disp('Forcing SDOnlyInPlace on for IQ recording');
          dopPRFNom = evParms.ev.dopPRFInPlace;
          UIValue = dopPRFNom;
          assignin('base', 'evParms', evParms);
          deferControlUpdate=0;
          ControlUpdate = changeprffn(UIValue, [], ...
                                      deferControlUpdate);
          
        end % evParms.state.SDOnlyInPlace | ~evParms.flag.forceHighPRFModeForRecordingIQ
      end % else of sdonlyinplace
      
      if 2==recordData.recordState
        if ~pausedMsgPrint
          disp([mfilename ': Recording paused']);
        end
        pausedMsgPrint=1;
      else
        pausedMsgPrint=0;                
      end     
    end % evParms.flag.forceRecord | (1 == successFlag)
    
    if switchToWT
      % this will run once all instances have stopped recording
      if ~ismember(evParms.state.recording, 1)
        % recording has stopped for all instances
        if evParms.flag.wtCapability
          evalin('base', 'evParms.flag.doWT = 1;');
        end
        if evParms.flag.BLineWTCapability
          evalin('base', 'evParms.flag.BLineWT = 1;');
        end
        
        evParms = evalin('base', 'evParms');         
        evParms = setoverlaysdpriorityfn(evParms);                  
        evParms.state.SDOnlyInPlace=0;
        disp('Forcing SDOnlyInPlace off for wall track mode');
        dopPRFNom = evParms.ev.dopPRF;
        UIValue = dopPRFNom;
        evParms.SD.enableInstances = [];
        assignin('base', 'evParms', evParms);
        deferControlUpdate=0;
        ControlUpdate = changeprffn(UIValue, [], ...
                                    deferControlUpdate);
        if isfield(evParms.UI, 'handleSDFig1') && ...
           ~isempty(evParms.UI.handleSDFig1) && ...
           ishandle(evParms.UI.handleSDFig1)
          delete(evParms.UI.handleSDFig1);
          evParms.UI.handleSDFig1 = [];
        end
        
        % close SD figures
        if isfield(evParms.UI, 'handleSDFig2') && ...
              ~isempty(evParms.UI.handleSDFig2) && ...
              ishandle(evParms.UI.handleSDFig2)
          delete(evParms.UI.handleSDFig2);
          evParms.UI.handleSDFig2 = [];
        end
        
      end
    end
    
    
  end % else of doInit
end % function saveIQWithTimeStamp
