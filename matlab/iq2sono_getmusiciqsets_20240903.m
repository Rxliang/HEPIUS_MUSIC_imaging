vtData = getenv('MUSIC_DATA');
% survival studies, 202408
%set_list_file = 'getmusiciqsets_20240903.mat';
%load(fullfile(vtData, set_list_file), 'music_sets');

map_file = 'offlinesdlarge_matfilemap_getmusiciqsets_20240904.mat';
load(fullfile(vtData, map_file), 'map');

% options for spectrogram
spec = [];
spec.HPFOrder=6;
spec.cutHPF_Hz=50;
spec.msWinLen = 50;
spec.winOffsetFrac = 0.05;
%spec.spectralMethod = 'burg';   spec.binThresh=8e-14;
spec.spectralMethod = 'fft';   spec.binThresh=9e-7;
spec.posPGramFlag = 0;
spec.isLog = 1;
spec.isGray = 1;

  
for ss = 1:length(map.matFile)
  if isempty(map.matFile{ss})
    continue
  end
  thisMatFile = map.matFile{ss};
  dset = load(thisMatFile);
  inPath = dset.inPath;
  inFilePrePre = dset.music_sets.set_names{dset.setInd}.name
  sgSel = dset.sgSel; % opportunity to modify segments used
  % opportunity to modify pixel indicies used is here:
  topInd = dset.topInd;
  iq2sono

  if ~isfield(map, 'iq2sono')
    map.iq2sono = [];
  end
  
  map.iq2sono{ss} = outMatFileBase;
  save(fullfile(vtData, map_file), 'map');
  
end

  
return

