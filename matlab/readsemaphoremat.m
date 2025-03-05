vtData = getenv('MUSIC_DATA');
if ~exist(vtData, 'dir')
  mkdir(vtData);
end

handles.tempDir = [vtData '/tmp/'];
if ~exist(handles.tempDir, 'dir')
    mkdir(handles.tempDir);
end

recordDataFileName = fullfile(handles.tempDir, 'recordDataFile.mat');
load(recordDataFileName)
recordData

