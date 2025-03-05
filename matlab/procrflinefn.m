function procrflinefn(RData) 

persistent myHandle handleWT cnt MMat colsM indVec filtParms fidRF ...
    rowInd numLinesPerFrame instanceNumber fileVersion 
persistent zAx_mm tAx_s COMP zSelInd tSelInd
persistent dateNumVec frameStore numFrameStore ...
           frameCount numOriginPerInstance
persistent vtData

global recordDataFileName; 
global recordData;
global semaphoreKey; semaphoreKey = 24537;
persistent semaphoreOpenedState;
persistent pausedMsgPrint;
persistent sessionName nframesWatchdog sessionNameRec stepNo stepTrialNo stepSegNo;

%extclockfn

saveRF=0; % old save method
pltRF=0;
pltM=1;
pltWT=0;

evParms = evalin('base', 'evParms');

% only accommodates single origin at present
originNo=1;
%size(RData)
%RDataCopy = RData;
RDataSel = RData(:, evParms.BLineWT.channels{originNo});

if ~evParms.flag.BLineWT
  % waiting for non-WT sequence to become active following switch
  % to SD
  return
end

if isempty(vtData)
  instanceNumber = 1;
  fileVersion=1;
  pausedMsgPrint=0;
  nframesWatchdog=0;
  vtData = getenv('MUSIC_DATA');
  if isempty(vtData)
    error('MUSIC_DATA not in environment');
  end
end

if isempty(recordDataFileName)
  recordDataPath = fullfile(vtData, 'nicp', 'tmp');
  if ~exist(recordDataPath, 'dir')
    error(['recordDataPath: ' recordDataPath ' does not exist!'])
  end  
  recordDataFileName = fullfile(recordDataPath, 'recordDataFile.mat');
end

if isempty(sessionName)
  sessionName = [];
end

if isempty(semaphoreOpenedState)
  [dum, osType] = unix('echo $OSTYPE');
  if dum==0 & (strcmp(osType(1:6), 'darwin') | strcmp(osType(1:5), ...
                                                      'linux'))
    nonWindowsHost=1;
    semaphoreOpenedState = 0;
  else
    nonWindowsHost=0;
    semaphore('create',semaphoreKey,1); %Initial semaphore value 1;
    semaphoreOpenedState = 1;
  end
end

successFlag = getrecordinginfofn(semaphoreOpenedState);

%recordData.stepType
%recordData.UIModePriority

if isfield(recordData, 'stepType') & ~isempty(recordData.stepType) & (recordData.stepType(1) == 'd') & ~recordData.UIModePriority
  evParms.flag.BLineWT = 0;
  finishRecording=1;
  exitWT = 1;
  disp('Exiting line wall track mode. Will write any open output file.')
else
  finishRecording=0;
  exitWT = 0;
end

Receive = evalin('base','Receive');
%rcvInd=257;
% be careful here:
rcvInd = evParms.BLineWT.TXIndBWT;

% these are in base wavelengths for Trans.Freq
endDepth_wvl = Receive(rcvInd(1)).endDepth;
keepDepth_wvl = min(evParms.BLineWT.MMaxDepth_mm/...
                    evParms.gate.lambda_mm, endDepth_wvl); 

nTimeSamp = Receive(rcvInd(1)).endSample-Receive(rcvInd(1)).startSample+ ...
    1;
nWvlRoundTrip = nTimeSamp/Receive(rcvInd(1)).samplesPerWave;
nWvl = nWvlRoundTrip/2;

maxRowsKeep = floor(keepDepth_wvl/nWvl*nTimeSamp);

%rowInd = Receive(rcvInd).startSample:Receive(rcvInd).endSample;

lumenLateral_mm = evParms.gate.x(instanceNumber)*evParms.gate.lambda_mm;
lumenDepth_mm = evParms.gate.z(instanceNumber)*evParms.gate.lambda_mm;

if ~exist('MMat', 'var') | isempty(MMat)
  dateNumVec=[];
  frameCount = 0;
  numOriginPerInstance = 1; % single origin per instance for now
  numFrameStore =  evParms.BLineWT.numFrameStore;
  
  numLinesPerFrame = length(rcvInd);
  for q = 1:numLinesPerFrame
    rowInd(q,:) = Receive(rcvInd(q)).startSample: ...
        Receive(rcvInd(q)).startSample+maxRowsKeep-1;
  end
  
  % each frame gives one "M-line" per ensemble per origin
  colsM = ceil(evParms.BLineWT.lineRate_Hz * evParms.BLineWT.MStoreDuration_s);
  rowsM = size(rowInd,2);
  MMat = zeros(rowsM, colsM, 'single');
  frameStore = zeros(rowsM, colsM*numOriginPerInstance, 'int16');  
  indVec = 1:colsM;
  evParms.Trans.frequency = evalin('base', 'Trans.frequency');

  actualEndDepth_mm = ((maxRowsKeep-1)/ ...
                       Receive(rcvInd(1)).samplesPerWave)*evParms.gate.lambda_mm/2;
  zSpacing_mm = (evParms.gate.lambda_mm/2)/Receive(rcvInd(1)).samplesPerWave;
  zAx_mm = 0:zSpacing_mm:actualEndDepth_mm;
  
  %  zAx_mm = 0:linspace(0, evParms.BLineWT.MMaxDepth_mm, ...
%                    rowsM);
  tAx_s = linspace(0, evParms.BLineWT.MStoreDuration_s, colsM);

  % set channels to average based on current (initial) gate
  % position
  Trans = evalin('base', 'Trans');
  evParms = findelementsrflinefn(evParms, Trans);
end

numEnsemble = length(evParms.BLineWT.TXIndBWTPre);
RF = [];

%evParms.BLineWT.elements{originNo}

% note that RData can change during keyboard, so let's nail it down.

for q = 1:numEnsemble
  RF(:,q) = mean(RDataSel(rowInd(q,:), :),2).';  
end

% check saving mode
if evParms.flag.forceRecord | (1 == successFlag)
  % intention is to begin a new file upon resume
  doRecording =  (1 == recordData.recordState || ...
                  3 == recordData.recordState); 
  finishRecording = ((0 == recordData.recordState || ...
                      2 == recordData.recordState)& ...
                      1 == ...
                     evParms.state.recordingWT(instanceNumber));

  if ~doRecording & evParms.state.recordingWT(instanceNumber)
    disp('RF line wall track recording state switched to OFF.');
  end
  
  if doRecording & ~evParms.state.recordingWT(instanceNumber)
    % just started recording
     disp('RF line Wall track recording state switched to ON.');
     MMat = zeros(size(MMat)); % force clearing of MMat
     indVec = 1:colsM;
  end
  
  evalin('base', ...
         ['evParms.state.recordingWT(' num2str(instanceNumber) ') = ' ...
          num2str(doRecording) ';']);
  
  % Put things in buffer only in started or resumed state
  if evParms.flag.forceRecord | doRecording | finishRecording    
    if ~strcmp(recordData.sessionName, sessionName) 
      disp([mfilename ': Detected new session: ' recordData.sessionName]);
      disp([mfilename ': Old session was: ' sessionName]);
      MMat = zeros(size(MMat)); % force clearing of MMat
      indVec = 1:colsM;
      sessionName =  recordData.sessionName;
      stepSegNo=1;
      finishRecording = 0; % don't finish any previous recording
    end
    
    dateNum = datenum(now);
    dateStr = datestr(dateNum, 'yyyymmddHHMMSSFFF'); 
    frameCount = frameCount + 1;
    storeInd = ((frameCount-1)*evParms.numBWTTxEnsemblePerFrame)+...
        (1:evParms.numBWTTxEnsemblePerFrame);
    frameStore(:, storeInd) = RF;
    dateNumVec(1, frameCount) = now;
    saveBuffer = (frameCount==numFrameStore) ||  finishRecording;        

    nframesWatchdog = nframesWatchdog+1;
    nframesWatchdog = rem(nframesWatchdog-1, ...
                          evParms.ev.watchdogInterval_frames)+1;
    watchdogUpdate = (nframesWatchdog == ...
                      evParms.ev.watchdogInterval_frames);
    if watchdogUpdate & instanceNumber==1
      vtData = getenv('MUSIC_DATA');
      watchdogFilename = fullfile(vtData, evParms.ev.watchdogFilename);
      fidw = fopen(watchdogFilename, 'w');
      if fidw
        fprintf(fidw, dateStr);
        fclose(fidw);
%        disp('Wrote to watchdog file');
      else
        disp(['Error: Watchdog file: ' watchdogFilename ...
              ' could not be written.']);
      end                      
    end
          
    if saveBuffer 
      vtData = getenv('MUSIC_DATA');     
      % check if we need another file to complete recording
      % this step 
      contFlag = (recordData.stepNo==stepNo) & ...
          (recordData.stepTrialNo == stepTrialNo) & ...
          strcmp(sessionNameRec, ...
                 recordData.sessionName);
      sessionNameRec = recordData.sessionName
      
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
                          '_bwt_s' num2str(recordData.stepNo) ... 
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
      
      stepNo = recordData.stepNo;
      stepTrialNo = recordData.stepTrialNo;

      % prevent overwrites
      fileFixStr = '';
      cntFix=0;
      cont=1;
      while cont
        rfFile = [outFilePre fileFixStr '.rf'];
          if exist(rfFile, 'file')
            disp(['*** Warning: ' rfFile ' exists!']);
            cntFix=cntFix+1;
            fileFixStr = ['_fix' num2str(cntFix)]; 
          else
            cont=0;
          end
      end
        
      if ~exist(outPath, 'dir')
        disp(['*** Output path does not exist. Making it to ' ...
              'avoid crash. ***']);
        mkdir(outPath);
      end

      fidRF = fopen(rfFile, 'w');
      if ~fidRF
        error(['Could not open file: ' rfFile]);
      else
        disp(['Opened file: ' rfFile]);
      end
     
      frameLastCol = storeInd(end);
      
      fwrite(fidRF, fileVersion, 'uint8');
      fwrite(fidRF, frameLastCol, 'uint16');
      fwrite(fidRF, numEnsemble, 'uint16');      
      fwrite(fidRF, numLinesPerFrame, 'uint16');      
      fwrite(fidRF, evParms.BLineWT.frequency_MHz, 'single');
      fwrite(fidRF, evParms.BLineWT.frameRate_Hz, 'single');
      fwrite(fidRF, lumenLateral_mm, 'single');
      fwrite(fidRF, lumenDepth_mm, 'single');     
      fwrite(fidRF, evParms.gate.xSVWidth_mm(instanceNumber), 'single');      
      fwrite(fidRF, evParms.gate.zSVWidth_mm(instanceNumber), 'single');
      fwrite(fidRF, length(rowInd), 'uint16');
      fwrite(fidRF, zAx_mm, 'single');     
      fwrite(fidRF, dateNumVec, 'double');                 
      fwrite(fidRF, frameStore(:,1:frameLastCol), 'int16');

      fclose(fidRF);
      
      frameCount = 0;
      dateNumVec = [];

      if finishRecording
        disp('Finished recording to open files.')
      end
      
      if instanceNumber == 1
        % save ultrasound image window
        saveas(...
            evalin('base','Resource.DisplayWindow(1).figureHandle'), ...
            [outFilePre '.fig']);
      end
	    
    end % saveBuffer
          
    if 2==recordData.recordState
      if ~pausedMsgPrint
        disp([mfilename ': Recording paused']);
      end
      pausedMsgPrint=1;
    else
      pausedMsgPrint=0;                
    end     
  end  % evParms.flag.forceRecord | doRecording | finishRecording  
end % evParms.flag.forceRecord | (1 == successFlag) checking
    % recording mode

if saveRF
   if ~exist('fidRF', 'var') | isempty(fidRF)
     dateStr = datestr(now, 'yyyymmdd_HHMM');
     vtData = getenv('MUSIC_DATA');
     filenameRF = fullfile(vtData, 'rf', [mfilename '_' dateStr ...
                         '.rf']);
     fidRF = fopen(filenameRF, 'w');
     if ~fidRF
       error(['Could not open file: ' filenameRF]);
     else
       disp(['Opened file: ' filenameRF]);
     end
     
     fwrite(fidRF, evParms.BLineWT.frequency_MHz, 'single');
     fwrite(fidRF, evParms.BLineWT.frameRate_Hz, 'single');
     fwrite(fidRF, lumenLateral_mm, 'single');
     fwrite(fidRF, lumenDepth_mm, 'single');     
     fwrite(fidRF, length(rowInd), 'uint16');
     fwrite(fidRF, zAx_mm, 'single');     
   end
end

if ~exist('filtParms', 'var') | isempty(filtParms)
  fCutHigh_Hz = 0.5;
  fsSig_Hz = evParms.BLineWT.lineRate_Hz;
  fCutHighNorm = fCutHigh_Hz/(fsSig_Hz/2); % Hz
  [filtParms.BH,filtParms.AH] = butter(4, fCutHighNorm, 'high');
end
  
if ~exist('cnt', 'var') | isempty(cnt)
  cnt=1;
end

Trans = evalin('base', 'Trans');
XLim = [1 maxRowsKeep];

if pltRF
  if isempty(myHandle)||~ishandle(myHandle)
    figure('name','Receive Signal','NumberTitle','off');
    myHandle = axes('XLim', XLim, 'YLim',[-2048  2048], ...
                    'NextPlot','replacechildren');
  end
end

if pltM % dispdiv('cnt', cnt, 1) 
  % to reinitialize image, close last figure handle
  if ~isfield(evParms.UI, 'MImageFigureHandle') | ...
      isempty(evParms.UI.MImageFigureHandle) | ...
     ~ishandle(evParms.UI.MImageFigureHandle)
    evParms.UI.MImageFigureHandle = figure;
    figureName = evParms.UI.MImageFigureTitle;

    set(evParms.UI.MImageFigureHandle, ...
        'Position', evParms.figPos.usM, ...
        'Tag',figureName );  
    
%    set(evParms.UI.MImageFigureHandle, 'MenuBar','none');
%    set(evParms.UI.MImageFigureHandle, 'Toolbar','none');
    assignin('base', 'evParms', evParms);
  end
  
  cMapLen = 256;
  % do reinitialize, close last image handle
  if ~isfield(evParms.UI, 'MImageImageHandle') | ...
      isempty(evParms.UI.MImageImageHandle) | ...
     ~ishandle(evParms.UI.MImageImageHandle)
    figure(evParms.UI.MImageFigureHandle);
    % z,t
    MSkip = [2 2];
    tSelInd = 1:MSkip(1):length(tAx_s);
    zSelInd = 1:MSkip(2):length(zAx_mm);
    evParms.UI.MImageImageHandle = image(tAx_s(tSelInd), zAx_mm(zSelInd), ...
                                         zeros(length(zSelInd), length(tSelInd)));

    %    CMap = parula(cMapLen);
    CMap = gray(cMapLen);
    colormap(CMap);
    xlabel('t (s)');
    ylabel('depth (mm)');
    
    assignin('base', 'evParms', evParms);
    freq_MHz = evParms.BLineWT.frequency_MHz;
    
    [atten, comp] = usattenfn(zAx_mm, freq_MHz, evParms.BLineWT.tissueAtten_dBpMHz);
    COMP = repmat(comp(:), 1, size(MMat,2));
  end
  
  % get the current (x,z) marker position
  zLumen_mm = evParms.gate.z(1)*evParms.gate.lambda_mm;
  xLumen_mm = evParms.gate.x(1)*evParms.gate.lambda_mm;
    
  if ~isfield(evParms.UI, 'MImageImageGateLineHandle') | ...
      isempty(evParms.UI.MImageImageGateLineHandle) | ...
     ~ishandle(evParms.UI.MImageImageGateLineHandle)
    zLine_mm = evParms.gate.zSVWidth_mm(instanceNumber)/2+zLumen_mm;
    evParms.UI.MImageImageGateLineHandle = line(gca(evParms.UI.MImageFigureHandle), ...
                                                tAx_s([1 end]), ...
                                                zLine_mm*[1 1], ...
                                                'linestyle', '--', ...
                                                'color', 'g');
  end
      
  for q = 1:numEnsemble
    MMat(:, indVec(1)) = RF(:,q);
    indVec = circshift(indVec,-1);
  end
    
%  [dum, indLat] = findclosestinvec(xAx_mm, xLumen_mm);
%  IQSel = IQIn(:, indLat);
 
%img = abs(MMat).* COMP;;
  img = log(sqrt(abs(MMat))).* COMP;;
  img = cMapLen*img/maxall(img);
  img = img(zSelInd, tSelInd);
  
  set(evParms.UI.MImageImageHandle, 'CData', img);
  
%  if dispdiv('cnt', cnt, 50)
%    size(get(evParms.UI.MImageImageHandle, 'CData'))
    %keyboard
%  end
    
  drawnow limitrate

end % pltM

% check if we are in WT mode
if exitWT
  evParms = setoverlaysdpriorityfn(evParms);                  
  evParms.state.SDOnlyInPlace=0;
  disp('Forcing SDOnlyInPlace off after wall track mode');
  dopPRFNom = evParms.ev.dopPRF;
  UIValue = dopPRFNom;
  evParms.SD.enableInstances = [1 2];
  assignin('base', 'evParms', evParms);
  deferControlUpdate=0;
  ControlUpdate = changeprffn(UIValue, [], ...
                              deferControlUpdate);
  if ~isempty(evParms.UI.MImageFigureHandle) & ...
      ishandle(evParms.UI.MImageFigureHandle)
    delete(evParms.UI.MImageFigureHandle);
    evParms.UI.MImageFigureHandle = [];
  end
  MMat=[]; % this will lead to reinit  
end

cnt=cnt+1;


end

