function [outFile, iSet] = readimfn(inFilePre, matInPath, fileVersion, ...
                                     overWrite)

if nargin < 4
  overWrite=1;
end

if nargin < 3
  fileVersion=1;
end

mfile = mfilename;
  
outFile = fullfile(matInPath, [mfile '_' inFilePre '.mat']);

if ~overWrite
  if exist(outFile, 'file') 
    %disp(['Output file: ' outFile ' exists!']);
    iSet = [];
    load(outFile, 'iSet');
    return
  end
end

if nargin < 2 | isempty(matInPath)
  vtData = getenv('MUSIC_DATA');
  matInPath = [vtData '/iqdata/'];
end

iSet = [];
iSet.matInPath = matInPath;

iSet.inFilePre = inFilePre;

fileRe = fullfile(matInPath, [inFilePre '.img']);
%matInPath
if ~exist(fileRe, 'file')
  error([fileRe ': does not exist']);
end

%fileRe
%lslrt(fileRe)

fid=[];
fid(1) = fopen(fileRe, 'rb');

head = [];
for i = 1:1
  head{i}.fileVersion = fread(fid(i), 1, 'uint64')';
  head{i}.sz = fread(fid(i), 3, 'uint64')';
  head{i}.prf = fread(fid(i), 1, 'double');
  head{i}.dopPRF = fread(fid(i), 1, 'double');
  head{i}.ttna_s = fread(fid(i), 1, 'double');
  head{i}.nDopPRIs = fread(fid(i), 1, 'double');
  head{i}.datenum = fread(fid(i), head{i}.sz(3), 'double');  
  head{i}.origin_wvl = fread(fid(i), 3, 'double');  
  head{i}.PDelta_wvl = fread(fid(i), 3, 'double');  
  head{i}.freq_MHz = fread(fid(i), 1, 'double');      
  head{i}.mode = char(fread(fid(i), 1, 'char'));
  head{i}.freqTW_MHz = fread(fid(i), 1, 'double');
end

fn = fieldnames(head{1});

for i = 1:length(fn)
  a = getfield(head{1}, fn{i});
  eval([fn{i} ' = a;'])
end

len = prod(sz);
nPRI = round(sz(3)/(ttna_s*dopPRF));

% maybe need also baseline offset

ttnaDop_s = ttna_s / (dopPRF*ttna_s);
framePeriod_s = nPRI*ttna_s;
iSet.c = 1540;

wvl2mm = iSet.c/head{1}.freq_MHz/1e3;
origin_mm = head{1}.origin_wvl*wvl2mm;
xDelta_mm = head{1}.PDelta_wvl(1)*wvl2mm;
yDelta_mm = head{1}.PDelta_wvl(2)*wvl2mm;
zDelta_mm = head{1}.PDelta_wvl(3)*wvl2mm;

% if we capture CData directly, we use PDelta(2) to store isotropic pixwid
if yDelta_mm == 0
  xAx_mm = origin_mm(1):xDelta_mm:origin_mm(1)+(xDelta_mm*(sz(2)-1));
  zAx_mm = origin_mm(3):zDelta_mm:origin_mm(3)+(zDelta_mm*(sz(1)-1));
else
  xAx_mm = origin_mm(1):yDelta_mm:origin_mm(1)+(yDelta_mm*(sz(2)-1));
  zAx_mm = origin_mm(3):yDelta_mm:origin_mm(3)+(yDelta_mm*(sz(1)-1));    
end

[X_mm, Z_mm] = meshgrid(xAx_mm, zAx_mm);

iSet.head = head{1};
iSet.nPRI=nPRI;
iSet.ttna_s = ttna_s;
iSet.wvl2mm = wvl2mm;
iSet.origin_mm = origin_mm;
iSet.xDelta_mm = xDelta_mm;
iSet.zDelta_mm = zDelta_mm;
iSet.X_mm = X_mm;
iSet.Z_mm = Z_mm;
iSet.xAx_mm = xAx_mm;
iSet.zAx_mm = zAx_mm;

if ~exist('im', 'var') | isempty(im)

  iqVec = fread(fid(1), len, 'single=>single');  
    
  if length(iqVec) < len
    error('Size mismatch');
  end
        
  im = reshape(iqVec, sz);

  iSet.im = im;

  save(outFile, 'iSet');
  lslrt(outFile)
end

end

