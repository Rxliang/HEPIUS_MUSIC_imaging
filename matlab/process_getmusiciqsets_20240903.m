vtData = getenv('MUSIC_DATA');

% survival studies, 202408
set_list_file = 'getmusiciqsets_20240904.mat';
load(fullfile(vtData, set_list_file), 'music_sets');
s = 31; % process this set in the list
inPath = fullfile(vtData, music_sets.set_names{s}.folder);
inFilePrePre = [music_sets.set_names{s}.name '_sg'];
inFilePre = [inFilePrePre]; 
segSet = 1:music_sets.set_names{s}.num_seg;
sgSel = segSet; % you can select a subset of time segments if needed

% if xzRect is not specified for this set:
%xzRectInitMode = 1; % search entire saved IQ region
xzRectInitMode = 2; % window around center of saved IQ region
xzRectWidth_mm = 0.5;
xzRectHeight_mm = 2;

% each set needs an xzRect in which to search for 
xzRect{1} = [-6 5
             6 8];

xzRect{2} = [0.37 7.7
             0.68 10 ];

xzRect{3} = [-0.80 7
             -0.65 8 ];

xzRect{4} = [0.4 4
             0.6 7];

xzRect{5} = [];

keepNTop = 10; % number of highest power pixels to store, and average in xzRect
               % the number actually used by iq2sono later on can be reduced

offlinesdlarge
