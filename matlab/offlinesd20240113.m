% copied fresh from gray vittamed

clear files spectralDoppler_jon_1

% need to have script vars in mem, Trans, Resource + ?
mfile = mfilename;

vtData = getenv('MUSIC_DATA');
vsExit = 0;

if 0
  inPath = [vtData '/iqdata/'];
%  inFilePrePre = 'sdc6s_20180823_0004_3'
  inFilePrePre = 'sdc6s_20180823_0004_3';
  inFilePre = ['iqim2sd_' inFilePrePre]; % ica
  
%  x_mm = -3;
%  z_mm = 64;

  x_mm = -1.7;
  z_mm = 67;
  
  xWid_mm=2;
  zWid_mm=2;
  
  xzOverride=1;
  
  imPlotInd = 500;
end


if 0
  inPath = [vtData '/iqdata/'];
%  inFilePrePre = 'sdc6s_20180822_2351_5';
  inFilePrePre = 'sdc6s_20180822_2351_4';
  inFilePre = ['iqim2sd_' inFilePrePre]; % ica

  x_mm = -5.4;
  z_mm = 56;
  
  xWid_mm=2;
  zWid_mm=2;
  
  xzOverride=1;
  
  imPlotInd = 200;
end


if 0
  inPath = [vtData '/iqdata/bcm20190305/'];
  inFilePrePre = 'bcm_test_20190226_20190305';
  inFilePre = ['iqim2sd_' inFilePrePre]; % ica

  x_mm = -5.4;
  z_mm = 56;
  
  xWid_mm=2;
  zWid_mm=2;
  
  xzOverride=1;
  
  imPlotInd = 200;
end

fileVersion=2;

mfileMat = 'readiqv2fn';

if 1
  inPath = [vtData '\nicp\'];
  inFilePrePre = 'subject1_20191025180720\subject1_s4_t1_sg1_i1';
  inPath = fullfile(vtData, 'nicp', 'semaphoretest');
  inFilePrePre = 'init_s1_t1_sg7_i1'
  inFilePre = [mfileMat '_' inFilePrePre]; % ica
  fileVersion=3;

  x_mm = -5.4;
  z_mm = 56;
  
  xWid_mm=2;
  zWid_mm=2;
  
  xzOverride=0;
  
  imPlotInd = 200;
end


matFile = [inPath inFilePre '.mat'];

if 1| ~exist(matFile, 'file')
  startTime_s=0;

  pwrThres = 0.035;
  DopState = 'computeCFIPowerEst';
  overWrite=1;
  matFile = readiqv2fn(inFilePrePre, startTime_s, inPath, fileVersion, ...
                       overWrite);
  lslrt(matFile)
  xzOverride=1;
  %inFilePre= ['readiqv2fn' inFilePrePre];
  Trans.frequency = 2.8409; 
end

load([fullfile(inPath, inFilePre) '.mat'], 'iqSet');

if xzOverride
  iqSet.x_mm = x_mm;
  iqSet.z_mm = z_mm;
  iqSet.xWid_mm=xWid_mm;
  iqSet.zWid_mm=zWid_mm;
end

xAx_mm = iqSet.X_mm(1,:);
zAx_mm = iqSet.Z_mm(:,1);

r = iqSet.head.rangeSubSampleIndex;
c = iqSet.head.lateralSubSampleIndex;

evParms.SD.dopWFCutoffNorm = 0.015; %allfilter normalized cutoff frequency
evParms.SD.SDDynRangeDB = 0.000013;% not including FFT processing gain
iqSet.evParms.SD = evParms.SD;
offlinesdfn20240113(iqSet, iqSet.IQMat, r, c);

  

