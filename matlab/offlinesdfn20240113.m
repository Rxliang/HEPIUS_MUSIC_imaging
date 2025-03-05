function offlinesdfn20240113(iqSet, IQMat, r, c);

evParms = [];
evParms.SD = iqSet.evParms.SD;

evParms.state.recording = 0;
evParms.state.SDOnlyInPlace = 1;
evParms.ev.nDopPRIs = iqSet.head.nDopPRIs;
evParms.gate.narrowbanding=0;
evParms.gate.rangeSubSampleIndex{1}=min(r):max(r);
evParms.gate.lateralSubSampleIndex{1}=min(c):max(c);
evParms.flag.offline=1;
evParms.ev.framePeriod_s = iqSet.nPRI*iqSet.ttna_s;
evParms.ev.dopPRF = iqSet.head.dopPRF;
evParms.flag.audioOn=0;
evParms.flag.instanceNoWithAudioEnabled=1;
evParms.state.CDIState=1; % irrelevant

evParms.ev.numSaveFramesPerBatch=1;

Trans.frequency = iqSet.head.freq_MHz;
nPulses = evParms.ev.nDopPRIs;
Resource.Parameters.speedOfSound = iqSet.c;
evParms.SD.titleFlow={'offline'};
evParms.SD.markerColorFlow = [ 1   1   0  ; 
                               0.2   1 .2];
evParms.SD.NWindow = round(iqSet.head.dopPRF*0.014)*2; 
evParms.SD.despeckle = 0.4;
%evParms.SD.dopWFCutoffNorm = 0.015*4; %allfilter normalized cutoff frequency
evParms.SD.dopSweepTimeIndex = 4%1-4 avail for 20msec frame time
evParms.SD.noisePersist=0.95;
%evParms.SD.SDDynRangeDB = 13;% not including FFT processing gain
%evParms.SD.SDDynRangeDB = 1;% not including FFT processing gain
%evParms.SD.SDDynRangeDB = 0.25;% not including FFT processing gain
%evParms.SD.SDDynRangeDB = 0.01;% not including FFT processing gain

evParms.SD.baseLineShiftNorm = -0.2;
evParms.SD.baseLineShiftNorm = -0.2;
evParms.SD.forceInit = 0;
evParms.flag.forceRecord=0;
assignin('base', 'evParms', evParms);
sgp= feval(@sdfn20240113,'spectrogramParameters',evParms.ev.framePeriod_s,...
           evParms.SD.dopSweepTimeIndex);
sweepTime = sgp.THist;
%baseLineShiftNorm = -0.375;

st=1;
en=nPulses;

%lastFrameNo = length(iqParms.datenumVec);

lastPulseNo = size(IQMat, 3)

while en <= lastPulseNo
  st
  en
%  spectralDoppler_jon_1(IQMat(:,:,st:en));
  % need to use function with any instanceNumber in suffix
  sdfn20240113(IQMat(:,:,st:en));
  pause(50e-3)
  st=en+1;
  en=en+nPulses;
%  pausede
end

  


