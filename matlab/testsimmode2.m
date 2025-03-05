vtData = getenv('MUSIC_DATA');
dataPath = fullfile(vtData);
%outFilePre = fullfile(dataPath, mfile);

% these are limited frame snapshots, not cine

fn = '20231214_MUSIC_Plane_Sag_WeakFlowVLow_Tubing.mat';
fn = '20231214_MUSIC_Plane_Axial_WeakVFlowLow.mat';
%fn = '20231214_MUSIC_FocusedDoppler_WeakFlowVLow.mat';
inMat = fullfile(dataPath, 'pig', ...
                 fn);

%if ~exist('d', 'var') | isempty(d)
%  d = load(inMat);
%end

%if ~exist('Receive', 'var') | isempty(Receive)
    %  d = load(inMat);
    %load(inMat, 'RcvData');
load(inMat);
              %end

Resource.Parameters.simulateMode=2;