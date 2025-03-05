vtData = getenv('MUSIC_DATA');
if ~exist(vtData, 'dir')
  mkdir(vtData);
end

handles.tempDir = [vtData '/tmp/'];
if ~exist(handles.tempDir, 'dir')
    mkdir(handles.tempDir);
end

global semaphoreKey; semaphoreKey = 24537;
global recordDataFileName; 
recordDataFileName = fullfile(handles.tempDir, 'recordDataFile.mat');
%Create semaphore - we use the mexfile associated with semaphore.c
semaphore('create',semaphoreKey,1); %Initial semaphore value 1;
recordData.subjectName = 'init'; 
dateStr = datestr(now, 'yyyymmdd_HHMM');
recordData.subjectName = dateStr;
recordData.recordState = 0;
recordData.sessionName = 'semaphoretest';
recordData.stepNo = 1;
recordData.stepTrialNo = 1;
recordData.stepType = ' '; % 'd'; %w'; %'d';
recordData.UIModePriority = 1;
%recordData.UIModePriority = 0; % this may set finishRecording 
%write data in shared file
semaphore('wait',semaphoreKey);
save(recordDataFileName,'recordData');
semaphore('post',semaphoreKey);

