function [outFile, iqSet] = readiqv2fn(inFilePre, startTime_s, ...
                                       matInPath, fileVersion, overWrite)

if nargin < 5
  overWrite=0;
end

if nargin < 2
  startTime_s=0;
end

if nargin < 4
  fileVersion=2;
end

mfile = mfilename;
  
outFile = fullfile(matInPath, [mfile '_' inFilePre '.mat']);

if ~overWrite & exist(outFile, 'file') 
  disp(['Output file: ' outFile ' exists!']);
  iqSet = [];
  return
end

if nargin < 3 | isempty(matInPath)
  vtData = getenv('VITTAMED_DATA');
  matInPath = [vtData '/iqdata/'];
end

iqSet = [];
iqSet.matInPath = matInPath;

iqSet.inFilePre = inFilePre;

fileRe = fullfile(matInPath, [inFilePre '.re']);
fileIm = fullfile(matInPath, [inFilePre '.im']);

%fileRe
%lslrt(fileRe)

if filesize(fileRe) ~= filesize(fileIm)
  error('Real and imag files not equal size');
end

fid=[];
fid(1) = fopen(fileRe, 'rb');
fid(2) = fopen(fileIm, 'rb');

head = [];
for i = 1:2
  head{i}.fileVersion = fread(fid(i), 1, 'uint64')';
  head{i}.sz = fread(fid(i), 4, 'uint64')';
  head{i}.prf = fread(fid(i), 1, 'double');
  head{i}.dopPRF = fread(fid(i), 1, 'double');
  head{i}.ttna_s = fread(fid(i), 1, 'double');
  head{i}.nDopPRIs = fread(fid(i), 1, 'double');
  head{i}.datenum = fread(fid(i), head{i}.sz(4), 'double');  
  if fileVersion > 1
    head{i}.origin_wvl = fread(fid(i), 3, 'double');  
    head{i}.PDelta_wvl = fread(fid(i), 3, 'double');  
    head{i}.freq_MHz = fread(fid(i), 1, 'double');      
  end

  if fileVersion > 2
    head{i}.mode = char(fread(fid(i), 1, 'char'));
    head{i}.freqTW_MHz = fread(fid(i), 1, 'double');
    head{i}.numLateral = fread(fid(i), 1, 'uint64');
    head{i}.lateralSubSampleIndex = fread(fid(i), head{i}.numLateral, ...
                                          'uint16');    
    head{i}.numRange = fread(fid(i), 1, 'uint64');
    head{i}.rangeSubSampleIndex = fread(fid(i), head{i}.numRange, 'uint16');        
  end
  

end

fn = fieldnames(head{1});

for i = 1:length(fn)
  a = getfield(head{1}, fn{i});
  b = getfield(head{2}, fn{i});  
  if a~=b
    error('Field mismatch between files');
  end
  eval([fn{i} ' = a;']);
end

len = prod(sz);
nPRI = round(sz(3)/(ttna_s*dopPRF));

% maybe need also baseline offset

ttnaDop_s = ttna_s / (dopPRF*ttna_s);

framePeriod_s = nPRI*ttna_s;
iqSet.c = 1540;

wvl2mm = iqSet.c/head{1}.freq_MHz/1e3;
origin_mm = head{1}.origin_wvl*wvl2mm;
xDelta_mm = head{1}.PDelta_wvl(1)*wvl2mm;
zDelta_mm = head{1}.PDelta_wvl(3)*wvl2mm;

xAx_mm = origin_mm(1):xDelta_mm:origin_mm(1)+(xDelta_mm*(sz(2)-1));
zAx_mm = origin_mm(3):zDelta_mm:origin_mm(3)+(zDelta_mm*(sz(1)-1));

[X_mm, Z_mm] = meshgrid(xAx_mm, zAx_mm);

iqSet.head = head{1};
iqSet.nPRI=nPRI;
iqSet.ttna_s = ttna_s;
iqSet.wvl2mm = wvl2mm;
iqSet.origin_mm = origin_mm;
iqSet.xDelta_mm = xDelta_mm;
iqSet.zDelta_mm = zDelta_mm;
iqSet.X_mm = X_mm;
iqSet.Z_mm = Z_mm;
iqSet.xAx_mm = xAx_mm;
iqSet.zAx_mm = zAx_mm;

if ~exist('iq', 'var') | isempty(iq)
  cnt=1;
  numFrames = 0;
  iq=[];

  while ~feof(fid(1)) 
    iqVec = fread(fid(1), len, 'single=>single');  
    
    if length(iqVec) < len
      break;
    end
    
    iqVecImag = fread(fid(2), len, 'single=>single');  
    if length(iqVecImag) < len
      break;
    end
    
    iqVec = iqVec+j*iqVecImag;
    
    szBlock = [sz(1) sz(2) sz(3)*sz(4)];
    
    tmp = reshape(iqVec, szBlock);
    
    blockLen = sz(3)*sz(4);
  
    iq(:,:, cnt:cnt+blockLen-1) = tmp;

    cnt = cnt+blockLen;
    
    numFrames = numFrames+sz(4); 

  end
end

if ndims(iq) == 4;
  IQMat = single(squeeze(iq)); % r+i*iqi);
else
  IQMat = single(iq);
end

szIQMat = size(IQMat);

lenIQ = szIQMat(end);

tAx = 0:ttnaDop_s:(lenIQ-1)*ttnaDop_s;
ind = find(tAx >= startTime_s);
IQMat = IQMat(:,:,ind);

iqSet.IQMat = IQMat;
iqSet.blockLen = blockLen;

%outFile = fullfile(matInPath, [inFilePre '.mat']);

if overWrite ~= -1
  save(outFile, 'iqSet', '-v7.3');
  lslrt(outFile)
end

end

