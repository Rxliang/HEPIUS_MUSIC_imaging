function saveimfn(imBuf)

persistent firstRun PDataThis
global recordData;

evParms = evalin('base', 'evParms');
semaphoreOpenedState = evParms.state.semaphoreOpenedState;
successFlag = getrecordinginfofn(semaphoreOpenedState);

% firstrun
if ~exist('firstRun', 'var') | isempty(firstRun)
  firstRun=0;
  PData = evalin('base', 'PData');                    
  %  PDataThis = PData(evParms.gate.PDataIndSDLarge);
  PDataThis = PData(evParms.gate.PDataIndB);
end

% if we have not had an active session and are not in an active recording state
% we should do nothing

% if any SD recording is ongoing, write images
doRecording = ~isempty(find(evParms.state.recording));

if  ~doRecording | ~successFlag | ...
    ( (0 == recordData.recordState || 2 == recordData.recordState) )
  return
end
  
% finish recording if we stop recording for any reason
% also, controling app should trigger a save of the past set
% by changing the step number

vtData = getenv('MUSIC_DATA');     
recordLocal = recordData;

outPath =  fullfile(vtData,  ...
                    recordLocal.sessionName);
            
outFilePrePre = fullfile(outPath, ...
                         [recordLocal.subjectName  ...
                    '_s' num2str(recordLocal.stepNo)]);            


if isfield(recordLocal, 'stepSegNo')
    outFilePre = [outFilePrePre ...
                  '_sg' num2str(recordLocal.stepSegNo)];
else
    outFilePre = outFilePrePre;
end

% if we already have an image of this step, return
fileExt = {'.img'};
imFile = [outFilePre fileExt{1}];

if exist(imFile, 'file')
    return
end
        
% prevent overwrites

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
imFile = [currentOutputFilePrefix '.img'];

instanceNumberUsed=1;
calledFn = 'saveim';
fileVersion=1;
                  
dopPRF = evParms.ev.dopPRF;
ttna = evParms.ev.ttna_microsec*1e-6;
PRF = evParms.ev.acqPRF;
freq = evParms.Trans.frequency;
closeFileAfterWriting=1;
            
feval(calledFn, single(imBuf), PRF, dopPRF, ttna, ...
      evParms.ev.nDopPRIs, ...
      imFile,  datenum(now), PDataThis.Origin, ...
      PDataThis.PDelta, freq, ...
      closeFileAfterWriting, instanceNumberUsed, ...
      char(evParms.this.mode), ...
      double(evParms.this.TWFreq_MHz), ...
      double(fileVersion));

%keyboard
end % function


