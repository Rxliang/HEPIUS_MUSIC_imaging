clear files spectralDoppler_jon_1 spectralDoppler

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

if 0
  inPath = [vtData '\nicp\iqdata\subject1_20191025180720'];
  inFilePrePre = 'subject1_s4_t1_sg1_i1';
  fileVersion = 2;
  %  inPath = [vtData '\nicp\semaphoretest\'];
  %inFilePrePre = 'init_s1_t1_sg7_i1'
  %inFilePre = [mfileMat '_' inFilePrePre]; % ica
  inFilePre = [inFilePrePre]; % ica
                              %  fileVersion=3;

  x_mm = -5.4;
  z_mm = 56;
  
  xWid_mm=2;
  zWid_mm=2;
  
  xzOverride=0;
  
  imPlotInd = 200;
end

if 0
  inPath = [vtData '/iqdata/'];
  inFilePrePre = 'bcm_test_20190226_20190305134843355_2';
  inFilePrePre = 'bcm_test_20190226_20190305134723847_2';
  
  x_mm = 0.5;
  z_mm = 48;
  %inFilePre = ['iqim2sd_' inFilePrePre]; % ica
  inFilePre = [inFilePrePre]; % ica
  fileVersion=2;
  %  x_mm = -5.4;
  %z_mm = 56;

  xWid_mm=2;
  zWid_mm=2;

  xzOverride=1;

  imPlotInd = 200;
  Trans.frequency = 2.8409;
end


if 0
  inPath = fullfile(vtData, 'iqdata', 'semaphoretest');
  
  inFilePrePre = '20240117_2033_s1_t1_sg1_i1';
  fileVersion = 3;
  %  inPath = [vtData '\nicp\semaphoretest\'];
  %inFilePrePre = 'init_s1_t1_sg7_i1'
  %inFilePre = [mfileMat '_' inFilePrePre]; % ica
  inFilePre = [inFilePrePre]; % ica
                              %  fileVersion=3;
  Resource.Parameters.speedOfSound = 1540;
  
  x_mm = -5.4;
  z_mm = 56;
  
  xWid_mm=2;
  zWid_mm=2;
  
  xzOverride=0;
  
  imPlotInd = 200;
end


if 1
    
  if 1
    inPath = fullfile(vtData, 'sss');
    inFilePrePre = 'sss20240307_003351_sdl_s3_sg1';
  end
        
    if 0
  inPath = [vtData '/pig/0215pig/0215Pigsurgery/'];
  inFilePrePre = '220240215_183529_s2_t1_sg4_i1'
   inFilePrePre = '220240215_183404_s1_t1_sg1_i1';
   end
   if 0
  inPath = [vtData '/pig/0215pig/02152024PigSugery/'];
  inFilePrePre = '120240215_182541_s1_t1_sg1_i1'
   end
   
  x_mm = 0.5;
  z_mm = 48;
  %inFilePre = ['iqim2sd_' inFilePrePre]; % ica
  inFilePre = [inFilePrePre]; % ica

  xWid_mm=2;
  zWid_mm=2;

  xzOverride=0;

  imPlotInd = 200;
  Trans.frequency = 2.8409;
end


matFile = [inPath inFilePre '.mat'];

if 1 | ~exist(matFile, 'file')
  startTime_s=0;
  pwrThres = 0.035;
  DopState = 'computeCFIPowerEst';
  overWrite=-1; % -1 is return iqSet, don't save to mat 
  fileVersion=3;
  matFile = readiqv2fn(inFilePrePre, startTime_s, inPath, fileVersion, ...
                       overWrite);
  lslrt(matFile)
  xzOverride=0;
  %inFilePre= ['readiqv2fn' inFilePrePre];
  %Trans.frequency = 2.8409; 
end

if 0
fileVersion = 1;
overWrite = 1;
[outFile, iSet] = readimfn(inFilePrePre, inPath, fileVersion, ...
                                     overWrite)

figure(10)
clf
for q = 1:size(iSet.im,3);
    imagesc(iSet.im(:,:,q).^0.25)
    colormap(gray)
    drawnow
end
end

load(matFile, 'iqSet');

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
  r = iqSet.head.rangeSubSampleIndex;
  c = iqSet.head.lateralSubSampleIndex;
end

%Trans.frequency = 2.8409
if fileVersion > 2
  Trans.frequency = iqSet.head.freq_MHz;
end


evParms.SD.dopWFCutoffNorm = 0.015; %allfilter normalized cutoff frequency
evParms.SD.SDDynRangeDB = 0.000013;% not including FFT processing gain
iqSet.evParms.SD = evParms.SD;
%offlinesdfn(iqSet, iqSet.IQMat, r, c);
offlinesd4p9p2fn(iqSet, iqSet.IQMat, r, c);

  

