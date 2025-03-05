mfile = mfilename;
addpath('icp-algorithm-3')

mfile = mfilename;

if ~isempty(time_sel_s{ss})
  startTime_s = time_sel_s{ss}(1);
  endTime_s = time_sel_s{ss}(2);  
else
  startTime_s=0;
  endTime_s = Inf;
end

overWrite=-1; % -1 is return iqSet, don't save to mat 
fileVersion=3;

% legacy manually selected sets for processing:
%setlist_iq2sono

if 0
  spec = [];
  spec.HPFOrder=6;
  spec.cutHPF_Hz=50;
  spec.msWinLen = 50;
  spec.winOffsetFrac = 0.05;
  spec.spectralMethod = 'burg';   spec.binThresh=8e-14;
  spec.spectralMethod = 'fft';   spec.binThresh=9e-7;
  spec.posPGramFlag = 0;
  spec.isLog = 1;
  spec.isGray = 1;
end


if isempty(dbstack('-completenames'))
  plt=1;
end

if length(topInd) == 1
  topIndStr = ['_ind' num2str(topInd)];
else
  topIndStr = ['_ind' vec2commadelim(topInd, 0, '_')];
end

outMatFileBase =  [mfile '_' spec.spectralMethod '_' inFilePrePre '_' vec2commadelim(sgSel, 0, '_') topIndStr];
outMatFilePre = fullfile(inPath, outMatFileBase);
outMatFile = [outMatFilePre '.mat'];

if exist(outMatFile, 'file')
  disp('Output file exists, skipping');
end

numSeg = length(sgSel);
inFilesPre = [];
for q = 1:numSeg
  inFilesPre{q} = [inFilePrePre '_sg' num2str(sgSel(q))];
end

sonoImFull = [];
timeAxFull_s = [];

iqSet = [];
lastTime_s = 0;
current_time_s = 0;

for q = 1:numSeg
  startTimeSeg_s = 0;
  [matFileDum, iqSets{q}] = readiqv2fn(inFilesPre{q}, startTimeSeg_s, inPath, ...
                                       fileVersion, overWrite);
  iqSet = iqSets{q};
  fs_Hz = iqSet.head.dopPRF;
  spec.T_s(1) = 1/fs_Hz; 
  sz = size(iqSet.IQMat);
  iqVecPre = reshape(iqSet.IQMat, [sz(1)*sz(2) sz(3)]);
  % mean across selected indices
  iqVec = [];
  %  iqVec{1} = double(mean(iqVecPre(380, :),1) - mean(iqVecPre(460, :),1));
  iqVec{1} = double(mean(iqVecPre(topInd, :),1));

  % wavelet analysis
  if 0
    len = length(iqVec{1});
    tAx_s = 0:1/fs_Hz:(len-1)/fs_Hz;
    [cfs,f] = cwt(iqVec{1}(:), fs_Hz);
 imagesc(timeAxSel_s, sonoSession.fAx,  log10(sonoImSel), im_lim);  
    hp = pcolor(tAx_s(:), log2(f), abs(cfs(:,:,2)));
    colormap(parula)
    hp.EdgeAlpha = 0;
    %  ylims = hp.Parent.YLim;
    %yticks = [-7:2]; %hp.Parent.YTick;
    %ylims = yticks([1 end]);
    cb = colorbar;
    cb.Label.String = 'Magnitude';
    axis tight
    hold on
    title('Magnitude scalogram')
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
    %  ylim(ylims);
 imagesc(timeAxSel_s, sonoSession.fAx,  log10(sonoImSel), im_lim);    %hp.Parent.YTick = yticks;
    %ytickformat('%1.1f')
    %hp.Parent.YTickLabels = 2.^yticks;
    hold off
  end
              
  [sonoSession, env, sumPowerEnv] = calcsonofn(iqVec, spec, plt);
  
  sig=1;
  img =  sonoSession.sonoWav.int;
  sonoPlot(sig).isLog = spec.isLog; 
  sonoPlot(sig).isGray = spec.isGray; 
  sonoPlot(sig).title = spec.spectralMethod; 
  sonoPlot(sig).img = img;
  sonoPlot(sig).env = [];

  if plt
  figure(sig);
  clf
  sonoSession.plotSono(sonoPlot);
  end

  %  prev_end_s = lastTime_s;
  %  if startTime_s >=  prev_end_s &  startTime_s >=  lastTime_s
  %    prev_end_s > endTime_s
  
  sonoImFull = cat(2, sonoImFull, single(sonoSession.sonoWav.int));
  timeAxFull_s = cat(2, timeAxFull_s, sonoSession.timeAx+lastTime_s);
  lastTime_s = lastTime_s + sonoSession.timeAx(end);  

  if lastTime_s >= endTime_s
    break
  end
  
end

ind_include = find( timeAxFull_s >= startTime_s & ...
                    timeAxFull_s <= endTime_s);

sonoImSel = sonoImFull(:,ind_include);
timeAxSel_s = timeAxFull_s(ind_include);

if plt
  figure(ff)
  nexttile
  %  imagesc(timeAxFull_s, sonoSession.fAx,  log10(sonoImFull), [10 20]);
  im_lim = [13 16];
  imagesc(timeAxSel_s, sonoSession.fAx,  log10(sonoImSel), im_lim);  
  colormap(gray)
  axis xy
  xlabel('time (s)')
  ylabel('frequency (Hz)')
  if analyze_nhood_flag
    ind_rc = ij2ind(r,c, sz(1));
    ind_str = num2str(ind_rc);
    title(['Gate at row ' num2str(r) ', column ' num2str(c) ', ind ' ind_str]);
  else
    ind_str = vec2commadelim(topInd);
    title(['Gate from indices: ' ind_str]);
  end
  

  title(tldc, replacechar(inFilePrePre, '_', '\_'))
end

spec.inFilesPre = inFilesPre;
spec.fs_Hz = fs_Hz;
spec.fAx_Hz = sonoSession.fAx;
spec.timeAxSel_s = timeAxSel_s;
spec.startTime_s = startTime_s;
spec.endTime_s = endTime_s;

%pausede
%save(outMatFile, 'sonoImFull', 'spec')
save(outMatFile, 'sonoImSel', 'spec')
lslrt(outMatFile)

