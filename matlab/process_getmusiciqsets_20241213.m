vtData = getenv('MUSIC_DATA');

% survival studies, 202408
set_list_file = 'getmusiciqsets_20241213.mat';
set_list_file = 'getmusiciqsets_20241227.mat';
load(fullfile(vtData, set_list_file), 'music_sets');

comment_m_file = 'music_dataset_comments_20241228';
eval(comment_m_file);


load(fullfile(vtData, set_list_file), 'music_sets');

%s = 31; % process this set in the list
s = 100; %first 1107 study, baseline
s = 101; %first 1107 study, closedloopfus
s=103; % 1107fus3, im and re files are 8192: wrong, why?
s=104;

s=1; % ok
s=2; % ok, short
s=3; % low signal, but something is there
s=4; % not useful, signal dies completely at end
s=5; % not useful
s=6; % not useful
s=7; % not useful, signal dies completely at end
s=8; % not useful, signal dies completely at end
s=9;  % poor, feedback pulsing, no cutoff
s=10; % poor, feedback pulsing, no cutoff
s=11; % poor, feedback pulsing, cuts off
s=12; % poor, feedback pulsing, cuts off
s=13; % very poor. may be somewhat redeemable. cuts off later on,
s=14; % very poor, cuts off 48 s. 
s=15; % useless
s=16; % looks like good signal ***
s=17; % looks like good signal ***
s=18; % strong signal, not at center of window
s=19; % strong signal, low velocity
s=20; % strong signal, low velocity
s=21; % strong signal, low velocity
s=22; % useless
s=23; % strong signal
s=24;% strong signal
s=25; % feedback, not good.
s=26; % useless
s=27; % poor, maybe not useless
s=28; %strong
s=29; % poor-fair
s=30; % poor, some adjust may help marginally, not done
s=31; % strong
s=32; % probably useless
s=33; % excellent
s=34; % useless
s=35; % useless
s=36; % useless
s=37; % some strong signal off center
s=38; % bad, feedback
s=39; % strong, but may not be consistent
s=40; % bad
s=41; % low, non-descript
s=42; % bad
s=43; % very strong signal, some aliasing
s=44; % useless
s=45; % useless
s=46; % useless
s=47; % poor. will convert anyway
s=48; % useless
s=49; % useless
s=50; % very strong, like 43
s=51; %useless
s=52; %useless
s=53; % some flat signal, but think useless
s=54; % some poor signal, converting anyway
s=55; % some poor flat signal, but think useless
s=56; % poor but something there, converting, not stable
s=57; % good, like 43,50
s=58; % also strong, but more feedback, middle more stable
s=59; % excellent
s=60; % has a mid-good signal
s=61; % has a some small area of signal, good B image, not very stable
s=62; % poor but process
s=63; % has a some small area of signal, good B image, not very stable
s=64; % B image weird, some strong const signals, but i feed useless
s=65; % B image weird, some strong const signals, but i feed useless
s=66; % useless
s=67; % useless, B image looks offset
s=68;  % looks uselesss non-pulsatile
s=69; % looks useless non-pulsatile
s=70; % strong signal, but not pulsatile
s=71; % moderate signal, but not pulsatile
s=72; % very good strong arterial
s=73; % has some pulsatile signal
s=74; % has some pulsatile signal. similar to 73
s=75; % has somewhat pulsatile signal
s=76; % has some weaker pulsatile signal
s=77; % not stable, keep but likely useless
s=78; % non-pulsatile, discard since likely useless
s=79; % non-pulsatile, discard since likely useless
s=80; % non-pulsatile, non-focal, discard since likely useless
s=81; % has some weaker pulsatile signal
s=82; % non-pulsatile, non-focal, discard since likely useless
s=83; % non-pulsatile, non-focal, discard since likely useless
s=84; % very good strong arterial
s=85; % very good strong arterial
s=86; % non-pulsatile, non-focal, discard since likely useless
s=87; % useless
s=88; % useless
s=89; % useless
s=90; % useless
s=91; % useless
s=92; % non-pulsatile, non-focal, discard since likely useless
s=93; % midling arterial
s=94; % strong arterial
s=95; % strong arterial
s=96; % strong arterial
s=97; % weak pulsatile (two areas), choosing stronger signal but less pulsatile
s=98; % good B, useless Dop
s=99; % useless
s=100; % 1107 good, long
s=101; % 1107 good
s=102; % 1107 good
s=103; % 103 - file trunc issue
s=104; % good
s=105; % good
s=106; % good
s=107; % good
       %s=108; % has some momentary pulsatile signal
       %s=109 % medium quality long recording. There is a more prox segment that is also viable
if 0
s=110; % first 1219 set useless
s=111; % useless
s=112; % error in sizes
s=113; % 1219_PostFusSession1_ClosedLoop, useless
s=114  % 1219_PostFusSession2_ClosedLoop, useless
s=115; % error in file
end

s=9;
s=26;
s=42;
s=32;
s=70;
s=71;
s=100;
%return
inPath = fullfile(vtData, music_sets.set_names{s}.folder);
inFilePrePre = [music_sets.set_names{s}.name '_sg'];
inFilePre = [inFilePrePre]; 
segSet = 1:music_sets.set_names{s}.num_seg;
sgSel = segSet; % you can select a subset of time segments if needed

sgSel=1:60;
% if xzRect is not specified for this set:
%xzRectInitMode = 1; % search entire saved IQ region
xzRectInitMode = 2; % window around center of saved IQ region

xzRectWidth_mm = 0.5;
xzRectHeight_mm = 2;

xzRect = [];

% each set needs an xzRect in which to search for 
xzRect{1} = [0.4 6.0
             0.6 7.5];

%xzRect{2} = [0.37 8.7
%             0.68 10 ];

%xzRect{3} = [-0.80 7
%             -0.65 8 ];

%xzRect{4} = [0.4 7
%             0.6 7.3];

xzRect{5} = [];

%xzRect{9} = [ 2.9 5.26
%              3.1 5.46];

%xzRect{9} = [ 2.92 5.44
%              2.93 5.45];

xzRect{9} = 'use_original';

xzRect{3} = [-0.9 7
             -0.6 7.4];

xzRect{13} = [2.0 5
              3.0 6.5];

xzRect{18} = [-2.2 6.5
              -1.9 7.0];

xzRect{21} = [2.2 8.0
              4   9.0];

xzRect{27} = [-1.78 6.4
              -1.30 7.2];

xzRect{31} = [-2.52 7.2
              -1.48 8.2];


xzRect{37} = [-2.4 7.0
              -1.85 7.32];

xzRect{43} = [2.5  8.85
              3.76 10.5];


xzRect{50} = [4.3  9.4
              4.62 10.2];


xzRect{60} = [-1.7  5.95
              -1.5  6.75];


xzRect{61} = [-1.77  5.72
              -1.72  6.0];


xzRect{63} = [-1.70  5.72
              -1.45  6.76];

xzRect{73} = [-1.82  6.02
              -1.53  6.50];


xzRect{75} = [1.85  8.23
              2.3   8.67];


xzRect{76} = [1.75  8.03
              2.2   8.78];


xzRect{77} = [1.9  7.98
              2.5  8.58];

xzRect{94} = [-3    8.026
              -2.1  8.913];

xzRect{95} = [3.37  9.31
              4.26  10.2];


xzRect{96} = [2.98  8.9
              4.00  10.25];

% weaker 97
%xzRect{97} = [-1.9  7.4
%              -1.6  7.7];

xzRect{97} = [-2.0   6.22
              -1.16  6.8];

xzRect{108} = [1.06   6.39
               1.36  7.42];


%xzRect{98} = [-0.97  6.96
%              -0.23  7.7];


keepNTop = 10; % number of highest power pixels to store, and average in xzRect
               % the number actually used by iq2sono later on can be reduced

overWriteIm = 0;
offlinesdlarge
