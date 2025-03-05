vtData = getenv('MUSIC_DATA');
comment = [];
ss=1; comment{ss} = 'ok';
ss=2; comment{ss} = 'ok, short';
ss=3; comment{ss} = 'low signal, but something is there';
ss=4; comment{ss} = 'not useful';
ss=5; comment{ss} = 'not useful';
ss=6; comment{ss} = 'not useful';
ss=7; comment{ss} = 'not useful';
ss=8; comment{ss} = 'not useful';
ss=9;  comment{ss} = 'poor, feedback pulsing, no cutoff';
ss=10; comment{ss} = 'poor, feedback pulsing, no cutoff';
ss=11; comment{ss} = 'poor, feedback pulsing, cuts off';
ss=12; comment{ss} = 'poor, feedback pulsing, cuts off';
ss=13; comment{ss} = 'very poor. may be somewhat redeemable.';
ss=14; comment{ss} = 'very poor, cuts off 48 s.';
ss=15; comment{ss} = 'useless';
ss=16; comment{ss} = 'looks like good signal';
ss=17; comment{ss} = 'looks like good signal';
ss=18; comment{ss} = 'strong signal, not at center of window';
ss=19; comment{ss} = 'strong signal, low velocity';
ss=20; comment{ss} = 'strong signal, low velocity';
ss=21; comment{ss} = 'strong signal, low velocity';
ss=22; comment{ss} = 'useless';
ss=23; comment{ss} = 'strong signal';
ss=24;comment{ss} = 'strong signal';
ss=25; comment{ss} = 'feedback, not good.';
ss=26; comment{ss} = 'useless';
ss=27; comment{ss} = 'poor, maybe not useless';
ss=28; comment{ss} = 'strong';
ss=29; comment{ss} = 'poor-fair';
ss=30; comment{ss} = 'poor, some adjust may help marginally, not done';
ss=31; comment{ss} = 'strong';
ss=32; comment{ss} = 'probably useless';
ss=33; comment{ss} = 'excellent';
ss=34; comment{ss} = 'useless';
ss=35; comment{ss} = 'useless';
ss=36; comment{ss} = 'useless';
ss=37; comment{ss} = 'some strong signal off center';
ss=38; comment{ss} = 'bad, feedback';
ss=39; comment{ss} = 'strong, but may not be consistent';
ss=40; comment{ss} = 'bad';
ss=41; comment{ss} = 'low, non-descript';
ss=42; comment{ss} = 'bad';
ss=43; comment{ss} = 'very strong signal, some aliasing';
ss=44; comment{ss} = 'useless';
ss=45; comment{ss} = 'useless';
ss=46; comment{ss} = 'useless';
ss=47; comment{ss} = 'poor. will convert anyway';
ss=48; comment{ss} = 'useless';
ss=49; comment{ss} = 'useless';
ss=50; comment{ss} = 'very strong, like 43';
ss=51; comment{ss} = 'useless';
ss=52; comment{ss} = 'useless';
ss=53; comment{ss} = 'some flat signal, but think useless';
ss=54; comment{ss} = 'some poor signal, converting anyway';
ss=55; comment{ss} = 'some poor flat signal, but think useless';
ss=56; comment{ss} = 'poor but something there, converting, not stable';
ss=57; comment{ss} = 'good, like 43,50';
ss=58; comment{ss} = 'also strong, but more feedback, middle more stable';
ss=59; comment{ss} = 'excellent';
ss=60; comment{ss} = 'has a mid-good signal';
ss=61; comment{ss} = 'has a some small area of signal, good B image, not very stable';
ss=62; comment{ss} = 'poor but process';
ss=63; comment{ss} = 'has a some small area of signal, good B image, not very stable';
ss=64; comment{ss} = 'B image weird, some strong const signals, but i feel it is useless';
ss=65; comment{ss} = 'B image weird, some strong const signals, but i feel is is useless';
ss=66; comment{ss} = 'useless';
ss=67; comment{ss} = 'useless, B image looks offset';
ss=68;  comment{ss} = 'looks uselesss non-pulsatile';
ss=69; comment{ss} = 'looks useless non-pulsatile';
ss=70; comment{ss} = 'strong signal, but not pulsatile';
ss=71; comment{ss} = 'moderate signal, but not pulsatile';
ss=72; comment{ss} = 'very good strong arterial';
ss=73; comment{ss} = 'has some pulsatile signal';
ss=74; comment{ss} = 'has some pulsatile signal. similar to 73';
ss=75; comment{ss} = 'has somewhat pulsatile signal';
ss=76; comment{ss} = 'has some weaker pulsatile signal';
ss=77; comment{ss} = 'not stable, keep but likely useless';
ss=78; comment{ss} = 'non-pulsatile, discard since likely useless';
ss=79; comment{ss} = 'non-pulsatile, discard since likely useless';
ss=80; comment{ss} = 'non-pulsatile, non-focal, discard since likely useless';
ss=81; comment{ss} = 'has some weaker pulsatile signal';
ss=82; comment{ss} = 'non-pulsatile, non-focal, discard since likely useless';
ss=83; comment{ss} = 'non-pulsatile, non-focal, discard since likely useless';
ss=84; comment{ss} = 'very good strong arterial';
ss=85; comment{ss} = 'very good strong arterial';
ss=86; comment{ss} = 'non-pulsatile, non-focal, discard since likely useless';
ss=87; comment{ss} = 'useless';
ss=88; comment{ss} = 'useless';
ss=89; comment{ss} = 'useless';
ss=90; comment{ss} = 'useless';
ss=91; comment{ss} = 'useless';
ss=92; comment{ss} = 'non-pulsatile, non-focal, discard since likely useless';
ss=93; comment{ss} = 'midling arterial';
ss=94; comment{ss} = 'strong arterial';
ss=95; comment{ss} = 'strong arterial';
ss=96; comment{ss} = 'strong arterial';
ss=97; comment{ss} = 'weak pulsatile (two areas), choosing stronger signal but less pulsatile';
ss=98; comment{ss} = 'good B, useless Dop';
ss=99; comment{ss} = 'useless';
ss=100; comment{ss} = '1107 good, long';
ss=101; comment{ss} = '1107 good';
ss=102; comment{ss} = '1107 good';
ss=103; comment{ss} = '103 - file trunc issue';
ss=104; comment{ss} = 'good';
ss=105; comment{ss} = 'good';
ss=106; comment{ss} = 'good';
ss=107; comment{ss} = 'good';
ss=108; comment{ss} = 'has some momentary pulsatile signal';
ss=109; comment{ss} = 'medium quality long recording. There is a more prox segment that is also viable';

map_file = 'offlinesdlarge_matfilemap_getmusiciqsets_20241227.mat';
load(fullfile(vtData, map_file), 'map');

map.comment = comment;


dateStr = datestr(now, 'yyyymmdd');
out_csv = fullfile(vtData, ['music_dataset_report_' dateStr '.csv']);
out_xlsx = fullfile(vtData, ['music_dataset_report_' dateStr '.xlsx']);
fid = fopen(out_csv, 'w');
fprintf(fid, 'Set index, Set name, Spectrogram file, Comment\n');

tab = [];
for q = 1:length(map.matFile)
  [dum, set_name, dum] = fileparts(map.matFile{q});
  if ismember(q, map.sel)
    sono_file = map.iq2sono{q};
  else
    sono_file = '';
  end
  %  fprintf(fid, '%d, %s,  %s, %s\n', q, set_name, sono_file, map.comment{q});
  tab.set_index(q,1) = q;
  tab.set_name{q,1} = set_name;
  tab.comment{q,1} = map.comment{q};
  tab.sono_file{q,1} = sono_file;
end

tabl = struct2table(tab);
delete(out_xlsx)
writetable(tabl, out_xlsx);


%fclose(fid)
lslrt(out_xlsx);


