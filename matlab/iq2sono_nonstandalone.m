% run before offline sd invoked in offlinesdlargecine

mfile = mfilename;
addpath('icp-algorithm-3')

fs_Hz = iqSet.head.dopPRF;

spec = [];
spec.HPFOrder=6;
spec.cutHPF_Hz=50;
spec.msWinLen = 50;
spec.winOffsetFrac = 0.05;
spec.spectralMethod = 'burg';   spec.binThresh=8e-14;
spec.spectralMethod = 'fft';   spec.binThresh=9e-7;
spec.posPGramFlag = 0;

spec.T_s(1) = 1/fs_Hz; 
plt=1;

sz = size(iqSet.IQMat);

iqVecPre = reshape(iqSet.IQMat, [sz(1)*sz(2) sz(3)]);
% mean across selected indices
iqVec = [];
iqVec{1} = mean(iqVecPre(topInd, :),1);

% wavelet analysis
if 0
  len = length(iqVec{1});
  tAx_s = 0:1/fs_Hz:(len-1)/fs_Hz;
  [cfs,f] = cwt(iqVec{1}(:), fs_Hz);
  
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
  %hp.Parent.YTick = yticks;
  %ytickformat('%1.1f')
  %hp.Parent.YTickLabels = 2.^yticks;
  hold off
end
              
[sonoSession, env, sumPowerEnv] = calcsonofn(iqVec, spec, plt)

sig=1;
img =  sonoSession.sonoWav.int;
sonoPlot(sig).isLog = 1; 
sonoPlot(sig).isGray = 1; 
sonoPlot(sig).title = spec.spectralMethod; 
sonoPlot(sig).img = img;
sonoPlot(sig).env = [];
 
%y3 = TVL1denoise(y3, 1.5, 100);
%y4a = medfilt1(y4,3);
%sonoPlot(1).img = imBin; %imgExt1;

    
figure(sig);
clf
sonoSession.plotSono(sonoPlot);

return

tAx = sonoSession.timeAx;
fAx = sonoSession.fAx;

figure(2)
clf
imagesc(tAx, fAx, img, [-2.16 6.8725])
%imagesc(log(img)) %, [-2.22 4.2659]);
%imagesc(img.^0.5, sqrt([-2.22 4.2659]))
colormap(gray)
axis xy
xlim([123 132])
ylim([0 1000])
xlabel('time (s)')
ylabel('frequency (Hz)');

return

imgCrop = img(1:50,50000:53000);
imagesc(imgCrop)
%envelope = sonoSession.getEnvelopeGMM(imgCrop);
%envelope = sonoSession.getEnvelope1(imgCrop);
%envelope = sonoSession.getSimpleEnv(imgCrop, 1e-3);
%envelope = sonoSession.envFilt(imgCrop);
%envelope = sonoSession.snsiEnvelope(imgCrop);
