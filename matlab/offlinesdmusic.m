%r = 1:4;
%c = 1:12;

r = 23:25;
c = 1:3 % 121:133; % pulsatile

c = 5:7;

iqSet = [];
evParms.SD = [];
evParms.SD.dopWFCutoffNorm = 0.015; %allfilter normalized cutoff frequency
evParms.SD.SDDynRangeDB = 0.000013;% not including FFT processing gain
                                   %iqSet.evParms.SD = evParms.SD;


iqSet.head.nDopPRIs =  size(IData{interIndSD},4);
% nPri is used for calculating framePeriod, so it should include all angles
numFrames = size(IData{interIndSD}, 5);
iqSet.nPRI = iqSet.head.nDopPRIs*dop.numAngs*numFrames;
iqSet.ttna_s = 1/dop.PRF/dop.numAngs;
iqSet.head.dopPRF = dop.PRF; 
iqSet.head.freq_MHz = Trans.frequency;
iqSet.c = P.c;

IQMatPre = squeeze(IData{interIndSD} + j * QData{interIndSD});

iqSet.IQMat = reshape(IQMatPre, [size(IQMatPre,1) size(IQMatPre,2) ...
                    size(IQMatPre,3)*size(IQMatPre,4)]);

vsExit = 0;
offlinesdfn(iqSet, iqSet.IQMat, r, c);

  

