% Script to get started up again after clearing workspace

bNew = 1;
if bNew
    a = sono;
    a.msWinLen = 50;
    winOffsetFrac = 0.05;
    a.loadAdcFile;
    a.gensonofwdrevfn(1, 'fft');
    a.gensonofwdrevfn(2, 'fft');
    imgInt1 = a.sonoWav.int;
    imgExt1 = a.sonoWav.ext;
    a.gensonofwdrevfn(1, 'burg');
    a.gensonofwdrevfn(2, 'burg');
%    a.statCompMap(1);
%    a.statCompMap(2);
%     a.binarizeSono(1);
%     a.binarizeSono(2);
end

img = a.sonoWav.int;

%img = a.KLD.int;
%img = a.binarySono.int;
% for plotting
sonoPlot(1).isLog = 1; 
sonoPlot(1).isGray = 1; 
sonoPlot(1).title = 'Low quality image'; 
sonoPlot(1).img = imgInt1;
sonoPlot(1).env = [];

%a.binarizeSono(a.sonoWav.int, .65);
% 
% sonoPlot(2).isLog = 0; 
% sonoPlot(2).isGray = 0; 
% sonoPlot(2).title = 'Internal binary'; 
%sonoPlot(2).img = a.binarySono;
% 
%a.edgeDetect(a.binarySono);
% 
% sonoPlot(3).isLog = 0; 
% sonoPlot(3).isGray = 0; 
% sonoPlot(3).title = 'Internal Edge'; 
% sonoPlot(3).img = a.binEdge;
% 
%a.getEnvelopeFromBinaryImg(a.binEdge, 1);
%y1 = smooth(a.binEnv);
% y2 = tvdip(a.binEnv, 2, 0);
% sonoPlot(2).env = y1;
% sonoPlot(3).env = y2;

%y3 = a.getEnvelopeGMM(a.sonoWav.int, 1.96);
%y3 = a.getSimpleEnv(img, 18);
%y3 = medfilt1(y3,3);
%y3 = a.snsiEnvelope(img);
%y3a = a.envFilt(img);
%y3 = a.getEnvelope1(img);
%y3 = a.getEnvelopeGMM(img);
%sonoPlot(1).env = y1;

% params.threshold = 1.5;  % 2000; .001
% params.numIters = 100;
% img1 = log10(a.simpleFilter(img, 3, params));
imBin = a.binarizeSono(img, .0005);
y3 = a.getEnvelopeFromBinaryImg(imBin, 0);
y3 = TVL1denoise(y3, 1.5, 100);
%y4a = medfilt1(y4,3);
%sonoPlot(1).img = imBin; %imgExt1;
sonoPlot(1).env = y3;

img = a.sonoWav.ext;
sonoPlot(2).isLog = 1; 
sonoPlot(2).isGray = 1; 
sonoPlot(2).title = 'High quality image'; 
sonoPlot(2).img = imgExt1;
%sonoPlot(2).env = y4;
%y4 = a.getSimpleEnv(img, 15);
%y4a = a.envFilt(img);
%y4 = a.snsiEnvelope(img);
%y4 = a.getEnvelope1(img);
%y4 = a.getEnvelopeGMM(img);
imBin = a.binarizeSono(img, .00005);
y4 = a.getEnvelopeFromBinaryImg(imBin, 0);
y4 = TVL1denoise(y4, 1.5, 100);
%sonoPlot(2).img = imBin;
sonoPlot(2).env = y4;

% y5 = a.getEnvelope1(img1);
% y5 = TVL1denoise(y5, 1.5, 100)*100;
% sonoPlot(2).env = y5;

%k = a.getKalmanEnv(img, y3, 500);
%sonoPlot(1).env = k(1,:);

figure;
a.plotSono(sonoPlot);
pause; 

envSono.int.fft = a.createEnvSono(imgInt1, y3);
envSono.int.AR = a.createEnvSono(a.sonoWav.int, y3);
envSono.ext.fft = a.createEnvSono(imgExt1, y4);
envSono.ext.AR = a.createEnvSono(a.sonoWav.ext, y4);

mom1.int.fft = a.spectralMoment(envSono.int.fft, y3, 1);
mom1.int.AR = a.spectralMoment(envSono.int.AR, y3, 1);
mom1.ext.fft = a.spectralMoment(envSono.ext.fft, y4, 1);
mom1.ext.AR = a.spectralMoment(envSono.ext.AR, y4, 1);

mom2.int.fft = a.spectralMoment(envSono.int.fft, y3, 2);
mom2.int.AR = a.spectralMoment(envSono.int.AR, y3, 2);
mom2.ext.fft = a.spectralMoment(envSono.ext.fft, y4, 2);
mom2.ext.AR = a.spectralMoment(envSono.ext.AR, y4, 2);

sonoPlot(1).img = imgInt1;
sonoPlot(1).env = mom1.int.fft;
sonoPlot(2).img = a.sonoWav.int;
sonoPlot(2).env = mom1.int.AR;


figure;
a.plotSono(sonoPlot);
% a.plotEnvelopes(y4a, y3a);
