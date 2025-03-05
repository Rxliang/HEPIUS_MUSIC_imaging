clear files spectralDoppler_4p9p2_saveiq
clear cdiprocoffline_nxtv1p0p4

% version for large SD file, does not save a mat file

mfile = mfilename;
vtData = getenv('MUSIC_DATA');
vsExit = 0;

mfileMat = 'readiqv2fn';

outMatFile = fullfile(inPath, [inFilePre vec2commadelim(sgSel, 0, '_') '.mat']);
outMatCDI =  fullfile(inPath, [inFilePre vec2commadelim(sgSel, 0, '_') '_cdi.mat']);

% file that records the outMatFile for each set in a data folder.
% it gets updated with each set processed
outMatFolderFile = fullfile(vtData, [mfile '_matfilemap_' set_list_file]);

if ~exist(outMatFolderFile, 'file')
    map = [];
    map.matFile{s} = outMatFile;
    map.inPath = inPath;
    map.music_sets = music_sets;
else
    load(outMatFolderFile);
    if length(map.matFile)>=s
      if ~isempty(map.matFile{s})
        disp('Warning: mat file mapping already exists for this dataset, ctrl-c to abort')
        pausede
      end
    end
    map.matFile{s} = outMatFile;
end

save(outMatFolderFile, 'map')
lslrt(outMatFolderFile);

if ~exist('keepTop', 'var') | isempty(keepTop)
  keepNTop = 10; % number of highest power pixels to average in xzRect
end

% manually select set
% setlist_offlinesdlarge.m

numSeg = length(sgSel);
inFilesPre = [];
for q = 1:numSeg
  inFilesPre{q} = [inFilePrePre num2str(sgSel(q))];
end

disp('Processing still images.');
isIm = [];
iSet = [];
for q = 1:numSeg
  imFile = fullfile(inPath, [inFilesPre{q} '.img']);

  if exist(imFile, 'file')
    fileVersion = 1;
    [outFile, iSet{q}] = readimfn(inFilesPre{q}, inPath, fileVersion, ...
                               overWriteIm);
  end  
  isIm(q)=1;
end

r = [];
c = [];
startTime_s=0;
pwrThres = 0.035;
DopState = 'computeCFIPowerEst';
overWrite=-1; % -1 is return iqSet, don't save to mat 
fileVersion=3;

if ~exist('iqSet', 'var') | isempty(iqSet)
  iqSet = [];

  for q = 1:numSeg
    [matFileDum, iqSets{q}] = readiqv2fn(inFilesPre{q}, startTime_s, inPath, ...
                                 fileVersion, overWrite);
    if q==1
      iqSet = iqSets{1};
      szi = size(iqSet.IQMat);
      tmp = iqSet.IQMat;
      disp('allocating memory for full IQMat');
      iqSet.IQMat = zeros(szi(1), szi(2), szi(3)*numSeg, 'single') + ...
          j*zeros(szi(1), szi(2), szi(3)*numSeg, 'single');
      disp('allocated memory for full IQMat');    
      iqSet.IQMat(:,:, 1:szi(3)) = tmp;
      blk_end = szi(3);
    else
      %    iqSet.IQMat = cat(3, iqSet.IQMat, iqSets{q}.IQMat);
      blk_start = (q-1)*szi(3)+1;
      blk_end = blk_start+size(iqSets{q}.IQMat,3)-1;
      iqSet.IQMat(:,:,blk_start:blk_end) = iqSets{q}.IQMat;
    end
    iqSet.set_start_datenum(q) = iqSets{q}.head.datenum(1);
  end
end % if iqSets exists

iqSet.inFilesPre = inFilesPre;
iqSet.blk_end = blk_end;

%for q = 1:numSeg
%  iqSet.set_start_datenum(q) = iqSets{q}.head.datenum(1);
%end

%iqSet.IQMat = iqSet.IQMat(:,:,1:blk_end);

iqSet.blockLen = size(iqSet.IQMat,3);
pwrThres = 0.15;
DopState = 'computeCFIPowerEst';
q=1;
np=20; % CDI pulses
ttnaDop_s =  1/iqSet.head.dopPRF;
frameDelay_s = np*ttnaDop_s;

if ~exist('imCDI', 'var') | isempty(imCDI)
  if exist(outMatCDI, 'file')
    disp(['Loading existing generated B and CDI images from: ' outMatCDI]);
    load(outMatCDI, 'imB', 'imCDI', 'pwr', 'mn', 'mx', 'tAx_s');
  else
    cnti=1;
    pwr = [];
    tAx_s = 0;
    mn = [];
    mx = [];
    imB = [];
    imCDI = [];
    while q+np-1 <= size(iqSet.IQMat,3)
      qStart = q;
      qEnd = qStart+np-1;
      iqm = iqSet.IQMat(:,:,qStart:qEnd);
      q=qEnd+1;
      sziqm = size(iqm);
      
      iqmThis = reshape(iqm, [sziqm(1:2) 1 np]); 
      
      imB(:,:,cnti) = single(abs(squeeze(mean(iqmThis,4))));  
      
      [~, imPre] = cdiprocessing_nxtv1p0p4(real(iqmThis), imag(iqmThis));
      imCDI(:,:,cnti) = single(abs(imPre));
      %im(:,:,cnti) = abs(cdiprocessing_jon(iqmThis));  
      
      pwr(cnti) = sumall(imCDI(:,:,cnti), 'all'); 
      mn(cnti) = minall(imCDI(:,:,cnti));
      mx(cnti) = maxall(imCDI(:,:,cnti));  
      
      cnti=cnti+1;
      tAx_s(cnti) = tAx_s(cnti-1)+frameDelay_s;    
    end % cdi process for np pulses
    disp('Saving B and CDI images');
    save(outMatCDI, 'imB', 'imCDI', 'pwr', 'mn', 'mx', 'tAx_s');
    lslrt(outMatCDI);
  end % load or regen cdi
end % imCDI exists


if 0
tAx_s = tAx_s(1:length(pwr));
figure(2)
clf
plot(tAx_s,pwr)

winLow = median(mn);
winHigh = median(mx);

if 1
figure(10)
clf
colormap gray
for q = 1:size(im,3)
  if 1
      imagesc(im(:,:,q), [winLow winHigh]);
      drawnow
      pause(frameDelay_s);
      pause
  end
end
end

end

im = imB;
numIm = size(im,3);
nDopPRIs = iqSet.head.nDopPRIs;
iqSampTime_s = 1/iqSet.head.dopPRF;

if 0
figure(10)
clf
if 0
skip = 100;
for fr = 1:skip:size(im,3)
  imagesc(im(:,:,fr))
  drawnow
  pause(iqSampTime_s*skip);
end
else
  imagesc(mean(im,3))
end
end

xzOverride=0;
if xzOverride
  iqSet.x_mm = x_mm;
  iqSet.z_mm = z_mm;
  iqSet.xWid_mm=xWid_mm;
  iqSet.zWid_mm=zWid_mm;
end

xAx_mm = iqSet.X_mm(1,:);
zAx_mm = iqSet.Z_mm(:,1);

if iqSet.head.fileVersion < 3
  [~,r] = findclosestinvec(zAx_mm, iqSet.z_mm);
  [~,c] = findclosestinvec(xAx_mm, iqSet.x_mm);
else
    %r = iqSet.head.rangeSubSampleIndex;
    %c = iqSet.head.lateralSubSampleIndex;
end
    
Trans.frequency = iqSet.head.freq_MHz;


%evParms.SD.dopWFCutoffNorm = 0.015; %allfilter normalized cutoff frequency
%evParms.SD.SDDynRangeDB = 0.000013;% not including FFT processing gain

evParms.SD.NWindow = round(iqSet.head.dopPRF*0.014)*2; 
evParms.SD.despeckle = 0.4;
evParms.SD.dopSweepTimeIndex = 1;%1-4 avail for 20msec frame time
evParms.SD.noisePersist=0.95;
%evParms.SD.SDDynRangeDB = 13;% not including FFT processing gain
evParms.SD.baseLineShiftNorm = 0;
evParms.SD.SDDynRangeDB = 6;
evParms.SD.dopWFCutoffNorm = 0.025;
iqSet.evParms.SD = evParms.SD;


figNoContext = 2;
figNoIm = 1;
stopFlag = 0;

showIm = mean(im,3);  
flowIm = mean(imCDI,3);
flowImNorm = flowIm/maxall(flowIm);
flowImRGB = repmat(flowImNorm, [1 1 3]);
flowImOvl = flowImRGB .* permute([1 0 0], [3 1 2]);

imRGB = repmat(showIm/maxall(showIm), [1 1 3]);

imBCDI = imRGB + 0.3*flowImOvl;

viewWinHalfWid_mm = 3;
viewWinHalfHt_mm = 3;
imBFull = [];
imCompNorm = [];
xViewCntr_mm = mean(iqSet.xAx_mm);
zViewCntr_mm = mean(iqSet.zAx_mm);

xlm = xViewCntr_mm + viewWinHalfWid_mm * [-1 1];
zlm = zViewCntr_mm + viewWinHalfHt_mm * [-1 1];
  
for q = 1:numSeg
  flowImNormInterp = interp2(iqSet.X_mm, iqSet.Z_mm, flowImNorm, iSet{q}.X_mm, ...
                             iSet{q}.Z_mm, ...
                             'bicubic', 0);
  flowImRGB = repmat(flowImNormInterp, [1 1 3]) .* permute([1 0 0], [3 1 2]);
  imBFull(:,:,q) = iSet{q}.im(:,:,1)/maxall(iSet{q}.im(:,:,1));
  imBFullRGB = repmat(imBFull(:,:,q), [1 1 3]);
  imCompNormPre = imBFullRGB+flowImRGB;
  imCompNorm(:,:,:,q) = imCompNormPre/maxall(imCompNormPre);
end

if 0
figure(figNoContext)
clf
tldc = tiledlayout(numSeg,2,'TileSpacing','Compact','Padding','Compact');
for q = 1:numSeg
  nexttile
  imagesc(iSet{q}.xAx_mm, iSet{q}.zAx_mm, imBFull(:,:,q));
  colormap gray
  nexttile
  imagesc(iSet{q}.xAx_mm, iSet{q}.zAx_mm, imCompNorm(:,:,:,q));
  xlabel('mm');
  ylabel('mm');
  ylim([0 20]);
end
end

% find xzRect_orig: original gate

if ~isfield(iqSet.evParms.SD, 'narrowbandInd')
  % grandfather in files that don't store this central pixel index
  szsd = size(iqSet.X_mm);
  ind_mat = reshape(1: prod(szsd), szsd);
  iqSet.evParms.SD.narrowbandInd = ind_mat(szsd(1)/2+1, szsd(2)/2);  
end

x_orig_mm = iqSet.X_mm(iqSet.evParms.SD.narrowbandInd);
z_orig_mm = iqSet.Z_mm(iqSet.evParms.SD.narrowbandInd);      
x_ind = find(iqSet.xAx_mm ==  x_orig_mm);
z_ind = find(iqSet.zAx_mm ==  z_orig_mm);

if strcmp(xzRect{s}, 'use_original')
  % use original gate           
  xo = x_orig_mm;
  zo = z_orig_mm;
  xl =  x_orig_mm-iqSet.xDelta_mm/2;
  xh =  x_orig_mm+iqSet.xDelta_mm/2;
  zl =  z_orig_mm-iqSet.zDelta_mm/2;
  zh =  z_orig_mm+iqSet.zDelta_mm/2;
  
  xzRect{s} = [xl zl;
               xh zh];

  xzRect{s} = [xo zo;
               xo zo];
end  

if ~exist('large_SD_dim', 'var')
  large_SD_dim = [];
end

large_SD_dim{s} = size(iqSet.X_mm);

% area restriction of search for highest power values
if ~exist('xzRect', 'var') | length(xzRect) < s | isempty(xzRect{s})
  % [xl zl
  %  xh zh]
  switch xzRectInitMode
    case 1       
      xl = iqSet.xAx_mm(1);
      xh = iqSet.xAx_mm(end);
      zl = iqSet.zAx_mm(1);
      zh = iqSet.zAx_mm(end);
    case 2
      xm = mean([iqSet.xAx_mm(1) iqSet.xAx_mm(end)]);
      zm = mean([iqSet.zAx_mm(1) iqSet.zAx_mm(end)]);            
      xl = xm-xzRectWidth_mm/2;
      xh = xm+xzRectWidth_mm/2;
      zl = zm-xzRectHeight_mm/2;
      zh = zm+xzRectHeight_mm/2;      
  end
        
  %  xzRect = [];
  xzRect{s} = [xl zl
               xh zh];
    
end

offlinesdlargecine

return

% use this to save the spectrogram. not needed if use iq2sono
sdOut.SDop.hFigure = [];
sdOut.SDop.HSpectrogram = [];

save(outMatFile, 'sdOut', 'topInd')
lslrt(outMatFile)
