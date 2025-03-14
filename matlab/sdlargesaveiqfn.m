function sdlargesaveiqfn(IQdataRe, IQdataIm)

persistent firstRun frameCount prevOutFilePrefix dateNumVec frameStore recordLocal PDataThis
global recordData;
global semaphoreKey;
global recordDataFileName;

% writing frames occurs when recording stops and there are frames in
% frameStore to write


if nargin == 2
  IQdata = IQdataRe + j*IQdataIm;
else
  IQdata = IQdataRe;
end

evParms = evalin('base', 'evParms');

if isfield(evParms, 'flag') &  isfield(evParms.flag, 'offline') & ...
        evParms.flag.offline ==1
  return
end

% maximum frames to store in memory before writing or losing frames.
numFramesStore = evParms.recording.numIQFrames;


semaphoreOpenedState = evParms.state.semaphoreOpenedState;
successFlag = getrecordinginfofn(semaphoreOpenedState);

% this needs to be evaluated each time because gate pos changes

PData = evalin('base', 'PData');                    
PDataThis = PData(evParms.gate.PDataIndSDLarge);

%keyboard
% firstrun
if ~exist('frameStore', 'var') | isempty(frameStore)  
  firstRun=1;
  sz = size(IQdata);
  switch length(sz)
    case 4
      frameStore = zeros(sz(1), sz(2), sz(4), numFramesStore, 'like', IQdata);
    case 2
      frameStore = zeros(sz(1), 1, sz(2), numFramesStore, 'like', IQdata);
    case 3
      frameStore = zeros(sz(1), sz(2), sz(3), numFramesStore, 'like', IQdata);
    otherwise
      error('Unaccommodated size of IQdata');
  end
      
  dateNumVec = zeros(1, numFramesStore, 'double');  
  
  % recording session init

  disp([mfilename ': doInit']);
  recordLocal.subjectName = recordData.subjectName;
  recordLocal.sessionName = recordData.sessionName;
  recordLocal.stepNo = recordData.stepNo;
  recordLocal.stepTrialNo = -1;    
  recordLocal.stepSegNo = 1;
  recordLocal.recordState = 0;
  % largeSD is always after last small SD gate
  recordLocal.instanceNumber = evParms.gate.numSD+1;  
  frameCount = 0;
end % firstrun

% if we have not had an active session and are not in an active recording state
% we should do nothing

if  ~successFlag | ...
    ( (0 == recordData.recordState || 2 == recordData.recordState) ...
      & (frameCount == 0) )
  return
end

% if we are not currently recording, we save whatever frames
% are in the store. otherwise we add frame to buffer until it is full

doRecording =  (1 == evParms.state.recording(recordLocal.instanceNumber)); 
  
% finish recording if we stop recording for any reason
% also, controling app should trigger a save of the past set
% by changing the step number

finishRecording=0;
sessionChange=0;
if doRecording 
  if ~strcmp(recordData.sessionName, recordLocal.sessionName) | ...
            (recordData.stepNo ~= recordLocal.stepNo) 
    % new session
    disp([mfilename ': Detected new session: ' recordData.sessionName ...
          ' step: ' num2str(recordData.stepNo)]);
    disp([mfilename ': Old session was: ' recordLocal.sessionName ...
          ' step: ' num2str(recordLocal.stepNo)]);
    sessionChange=1;
    if frameCount > 0 
      finishRecording=1; % save whatever is in buffer
    end
  end
  
  saveBuffer = (frameCount==numFramesStore) ||  finishRecording;
        
  % save buffer after set number of frames or when recording should end
  if saveBuffer 
    vtData = getenv('MUSIC_DATA');     
            
    outPath =  fullfile(vtData,  ...
                        recordLocal.sessionName);
            
    outFilePrePre = fullfile(outPath, ...
                             [recordLocal.subjectName '_sdl' ...
                              '_s' num2str(recordLocal.stepNo)]);            
      
    outFilePre = [outFilePrePre ...
                  '_sg' num2str(recordLocal.stepSegNo)];
      
    recordLocal.stepSegNo = recordLocal.stepSegNo+1;
      
    % prevent overwrites
    fileExt = {'.re', '.im'};
    fileFixStr = stopfileoverwritefn(outFilePre, fileExt);
                                                      
    if ~exist(outPath, 'dir')
      disp(['*** Output path does not exist. Making it to ' ...
            'avoid crash. ***']);
      mkdir(outPath);
    end
            
    evParms.this.mode = 'd';
    evParms.this.TWFreq_MHz = ...
        evParms.TW(evParms.ind.TWIndSD).Parameters(1);

    currentOutputFilePrefix = fileFixStr;
    iqFileRe = [currentOutputFilePrefix '.re'];
    iqFileIm = [currentOutputFilePrefix '.im'];

    instanceNumberUsed=1;
    %    if isreal(frameStore)
    %     keyboard
        %end    
    saveImFlag=1;
    saveiqfn(frameStore(:,:,:, 1:frameCount), evParms, PDataThis, ...
             iqFileRe, ...
             iqFileIm, dateNumVec(1,1:frameCount), ...
             instanceNumberUsed, saveImFlag);
    end % saveBuffer          
end % do or finish recording  

% if there is now a new session, update these
if sessionChange
  recordLocal.subjectName =  recordData.subjectName;
  recordLocal.sessionName =  recordData.sessionName;
  recordLocal.stepNo =  recordData.stepNo;
  recordLocal.stepSegNo=1;
  frameCount = 0;
else
  % update seg no
  if isfield(recordData, 'stepSegNo') & recordData.stepSegNo ~= recordLocal.stepSegNo;
    recordData.stepSegNo = recordLocal.stepSegNo;
    semaphore('wait',semaphoreKey);
    save(recordDataFileName,'recordData');
    recordDataFileName
    semaphore('post',semaphoreKey);
  end
  
end

% only add to buffer if we are recording
if (recordData.recordState==0 | recordData.recordState==2)
  return
end

frameCount = frameCount+1;
frameCount = rem(frameCount-1,numFramesStore)+1;

if frameCount==round(frameCount/10)*10
    disp(['Large SD frame ' num2str(frameCount) ' of ' ...
          num2str(numFramesStore)]);
end

dateNum = datenum(now);    
dateNumVec(1,frameCount) = dateNum;

sz = size(IQdata);
if ndims(IQdata) ==2
  iqSto = reshape(IQdata, sz(1), 1, sz(2));
else
  iqSto = squeeze(IQdata);
end
%keyboard
frameStore(:, :, :, frameCount) = iqSto;

firstRun=0;

end % function


