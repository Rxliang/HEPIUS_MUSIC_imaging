vtData = getenv('MUSIC_DATA');
music_home = getenv('MUSIC_HOME');
% survival studies, 202408
%set_list_file = 'getmusiciqsets_20240903.mat';
%load(fullfile(vtData, set_list_file), 'music_sets');
addpath(fullfile(music_home, 'matlab','icp-algorithm-3'));

screen_size = 1e3*[ 0.0010    0.0010    2.9257    1.2343];

%map_file = 'offlinesdlarge_matfilemap_getmusiciqsets_20241213.mat';
map_file = 'offlinesdlarge_matfilemap_getmusiciqsets_20241227.mat';
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
plt=1;
%return
sel = [1 2 3 9:14 16:21 23 24 28 29 30 31 33 37 39 43 47 50 54 57 58:63 70:76 77 81 84 85 93 94 95 96 97 100 101 102 104 105 106 107 108 109]
sel=107;
sel=9;
sel=26;
sel=42;
sel=32;
sel = 71;

sel=100;

time_sel_s = [];
time_sel_s{9} = [0 Inf];
time_sel_s{26} = [0 Inf];
time_sel_s{32} = [0 Inf];
time_sel_s{42} = [0 Inf];
time_sel_s{70} = [0 Inf];
time_sel_s{71} = [0 Inf];
time_sel_s{100} = [0 Inf];

nnhood = [];
nnhood{9} = 0; % how many nearest neighbors to process
nnhood{26} = 0; % how many nearest neighbors to process
nnhood{32} = 0; % how many nearest neighbors to process
nnhood{42} = 0; % how many nearest neighbors to process
nnhood{70} = 0; % how many nearest neighbors to process
nnhood{71} = 1; % how many nearest neighbors to process
nnhood{100} = 0; % how many nearest neighbors to process

topIndOverride = [];
topIndOverride{9} = [380 420 381 421];

topIndOverride{26} = 381;
topIndOverride{26} = [340 380]; % 340 better than 380 as center
topIndOverride{26} = [341 381]; % 341 better than 381 as center

topIndOverride{32} = 381;
topIndOverride{32} = 380;
topIndOverride{32} = [379 380];

topIndOverride{42} = [381];
topIndOverride{42} = [340];
topIndOverride{42} = [339 379 340 380];
%topIndOverride{42} = [379];

topIndOverride{70} = [381];
topIndOverride{70} = [301 341];
topIndOverride{70} = [341];

topIndOverride{71} = [381];

topIndOverride{100} = [];

for ss = sel % 104 % 1:length(map.matFile)
  if isempty(map.matFile{ss})
    continue
  end
  thisMatFile = map.matFile{ss};
  dset = load(thisMatFile);

  if ~isempty(topIndOverride{sel})
      topInd = topIndOverride{sel};
    else
      topInd = dset.topInd;
  end
  
  % if we have a single point gate, we can analyze its neighborhood
  analyze_nhood_flag = 0;
  if length(topInd)==1 & ~isempty(nnhood{sel}) & nnhood{sel} > 0
    szsd = dset.large_SD_dim{ss};
    ind_mat = reshape(1:prod(szsd), szsd);
    [r0,c0] = ind2ij(topInd, szsd(1));
    %  r0 = 4; c0=4;
    %  r0 = 4; c0=7;
    r_span = r0-nnhood{sel}:r0+nnhood{sel};
    c_span = c0-nnhood{sel}:c0+nnhood{sel};
    analyze_nhood_flag = 1;
  else
    r_span=1; % for figure tile layout
    c_span=1; % for figure tile layout    
  end

  inPath = dset.inPath;
  inFilePrePre = dset.music_sets.set_names{dset.setInd}.name
  sgSel = dset.sgSel; % opportunity to modify segments used
  % opportunity to modify pixel indicies used is here:
  ff=2; % figure index
  figure(ff)
  switch length(r_span)
    case 1
      pos = 1e3*[ 0.0153    0.5793    1.0651    0.5709];    
    case 3
      pos = 1e3*[ 0.0153    0.5793    1.0651    0.5709];
    case 5
      pos = 1e3*[ 0.0153    0.4724    1.4708    0.6778];
    case 7
      pos = 1e3*[ 0.0153    0.1250    1.8000    1.0252];
  end
  
  pos_rel_fig = pos ./ screen_size;
  pos_abs = setfigposrelfn(gcf, pos_rel_fig);

  tile_side = [length(r_span) length(c_span)];
  tldc = tiledlayout(tile_side(1),tile_side(2),'TileSpacing','Compact','Padding','Compact'); 

  if analyze_nhood_flag
    cnt = 0;
    for r = r_span
      for c = c_span
        topInd = ij2ind(r,c, szsd(1));
        iq2sono
        cnt=cnt+1;
        if r==r0 & c==c0
          outMatFileBase0 =  outMatFileBase;
        end
        map.iq2sono_nn{ss}{cnt} = outMatFileBase;  
      end
    end
  else    
    iq2sono
    outMatFileBase0 =  outMatFileBase;
  end
    
  if ~isfield(map, 'iq2sono')
    map.iq2sono = [];
  end
  map.sel = sel; 
  map.iq2sono{ss} = outMatFileBase0;
  save(fullfile(vtData, map_file), 'map');
  %  pausede
end

  
return

