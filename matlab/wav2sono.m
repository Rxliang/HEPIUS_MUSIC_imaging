mfile = mfilename;
addpath('icp-algorithm-3')

filename = '../data/filtered_240307.wav';
[aud44,fs44] = audioread(filename);

figure(1)
plot(aud44)

decFactor=10;
aud = decimate(aud44, decFactor);
fs = fs44/decFactor;

spec = [];
spec.HPFOrder=6;
spec.cutHPF_Hz=50;
spec.msWinLen = 50;
spec.winOffsetFrac = 0.05;
spec.spectralMethod = 'burg';   spec.binThresh=8e-14;
spec.spectralMethod = 'fft';   spec.binThresh=9e-7;

spec.T_s(1) = 1/fs; 
plt=1;

iqVecInt = aud(1:end);

iqVec = [];
iqVec{1}  = 0.5*(iqVecInt(:)+j*iqVecInt(:));

[sonoSession, env, sumPowerEnv] = calcsonofn(iqVec, spec, plt)

sig=1;
img =  sonoSession.sonoWav.int;
sonoPlot(sig).isLog = 1; 
sonoPlot(sig).isGray = 1; 
sonoPlot(sig).title = 'Low quality image'; 
sonoPlot(sig).img = img;
sonoPlot(sig).env = [];
 
%y3 = TVL1denoise(y3, 1.5, 100);
%y4a = medfilt1(y4,3);
%sonoPlot(1).img = imBin; %imgExt1;


    
figure(sig);
sonoSession.plotSono(sonoPlot);
    
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

imgCrop = img(1:50,50000:53000);
imagesc(imgCrop)
%envelope = sonoSession.getEnvelopeGMM(imgCrop);
%envelope = sonoSession.getEnvelope1(imgCrop);
%envelope = sonoSession.getSimpleEnv(imgCrop, 1e-3);
%envelope = sonoSession.envFilt(imgCrop);
%envelope = sonoSession.snsiEnvelope(imgCrop);