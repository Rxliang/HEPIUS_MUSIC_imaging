function [dataFileName]=saveimwithtimestampfn(x, PDataThis, ...
                                              NframesIn, instanceNumber,...
                                              semaphoreOpenedState)
% save image buffer frames
% intended to work in parallel with recording via sdfn.m
% it reads the same semaphore, but does not try to
% do any mode switching.
    
evParms = evalin('base', 'evParms');
verbose = 1;
dataFileName = '';
persistent firstRun X Timestamp TimestampNum nframes Nframes nframesWatchdog ...
    sessionName sessionNameRec pausedMsgPrint stepNo stepTrialNo stepSegNo;
  
persistent recordState

if isequal(x,'cleanup')
  [X, nframes, Nframes]=deal([]);
  disp('saveimwithtimestampfn: cleanup');
  nframesWatchdog=[];
  sessionName = '';
  sessionNameRec = ''; % tracked in saveBuffer conditional
  return
end

%Check for initialization

if ~exist('firstRun', 'var') | isempty(firstRun) 
  firstRun = 1;
end

if firstRun
    recordState = 0;
    disp('saveimwithtimestampfn: doInit');
    pausedMsgPrint=0;
    %initialization
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
    firstRun=0;
end

global recordData;
successFlag = getrecordinginfofn(semaphoreOpenedState);

if evParms.flag.forceRecord | (1 == successFlag)
  % intention is to begin a new file upon resume
  doRecording =  (1 == recordData.recordState || ...
                  3 == recordData.recordState); % start/resume
  
  % finish recording if we stop recording for any reason
  finishRecording = ((0 == recordData.recordState || ...
                      2 == recordData.recordState) & ...
                      1 == evParms.state.recording(instanceNumber));
          
  % Put things in buffer only in started or resumed state
  if evParms.state.SDOnlyInPlace | ...
     ~evParms.flag.forceHighPRFModeForRecordingIQ             
      if evParms.flag.forceRecord | doRecording | finishRecording    
        if ~strcmp(recordData.sessionName, sessionName) | ...
                (recordData.stepNo ~= stepNo)
          % new
          % session
          disp([mfilename ': Detected new session: ' recordData.sessionName ...
                ' step: ' num2str(recordData.stepNo)]);
          disp([mfilename ': Old session was: ' sessionName ...
                ' step: ' num2str(stepNo)]);
          nframes=0;
          X = zeros([size(x) Nframes], 'single');
          sessionName =  recordData.sessionName;
          stepNo =  recordData.stepNo;
          stepSegNo=1;
          finishRecording = 0; % don't finish any previous recording
        end

        nframes = nframes+1;
        nframes = rem(nframes-1,Nframes)+1;
        if nframes==round(nframes/20)*20
          disp(['Image frame ' num2str(nframes) ' of ' num2str(Nframes)]);
        end
        
        X(:,:,nframes) = single(x);
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
          vtData = getenv('MUSIC_DATA');
          watchdogFilename = fullfile(vtData, evParms.ev.watchdogFilenameIm);
          fidw = fopen(watchdogFilename, 'w');
          if fidw
              fprintf(fidw, dateStr);
              fclose(fidw);
          else
              disp(['Error: Watchdog file: ' watchdogFilename ...
                    ' could not be written.']);
          end                      
        end
        
        % save buffer after set number of frames or when recording should end
        if saveBuffer 
            vtData = getenv('MUSIC_DATA');     
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
            
            outPath =  fullfile(vtData,  ...
                                recordData.sessionName);
            
            outFilePrePre = fullfile(outPath, ...
                                  [recordData.subjectName ...
                                '_s' num2str(recordData.stepNo) ... 
                                '_t' ...
                                num2str(recordData.stepTrialNo)]);
            
            if instanceNumber == 1 & evParms.largeSDParms.largeSDSave
              % make available for SDLarge recording
              evalin('base', ['currentOutputFilePrefixIm = ''' outFilePrePre ...
                              ''';']);
            end
            
            outFilePre = [outFilePrePre ...
                          '_sg' num2str(stepSegNo) ...
                          '_i' num2str(instanceNumber)];  
            
            
            
            % prevent overwrites
            fileExt = {'.img'};
            fileFixStr = stopfileoverwritefn(outFilePre, fileExt);
            
            stepNo = recordData.stepNo;
            stepTrialNo = recordData.stepTrialNo;
                                                      
            if ~exist(outPath, 'dir')
              disp(['*** Output path does not exist. Making it to ' ...
                    'avoid crash. ***']);
              mkdir(outPath);
            end
            
            evParms.this.mode = 'd';
            evParms.this.TWFreq_MHz = ...
                evParms.TW(evParms.ind.TWIndSD).Parameters(1);

            %            nframes
            %TimestampNum(1:nframes)
            %size(X)
            %            keyboard
                               
            calledFn = 'saveim';
            fileVersion=1;
                  
            dopPRF = evParms.ev.dopPRF;
            ttna = evParms.ev.ttna_microsec*1e-6;
            PRF = evParms.ev.acqPRF;
            freq = evParms.Trans.frequency;
            closeFileAfterWriting=1;
            
            %            imFile
            feval(calledFn, X(:,:,1:nframes), PRF, dopPRF, ttna, ...
                  evParms.ev.nDopPRIs, ...
                  imFile, TimestampNum(1:nframes), PDataThis.Origin, ...
                  PDataThis.PDelta, freq, ...
                  closeFileAfterWriting, instanceNumber, ...
                  char(evParms.this.mode), ...
                  double(evParms.this.TWFreq_MHz), ...
                  double(fileVersion));


            % reset buffer
            X = zeros([size(x) Nframes], 'single');
            nframes=0;

	    
          end % saveBuffer
          
      end % evParms.flag.forceRecord | doRecording
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
    
end % function
