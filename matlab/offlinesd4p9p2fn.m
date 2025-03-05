function [sgp] = offlinesd4p9p2fn(iqSet, IQMat, r, c, init)

persistent SDop st 

if nargin < 5
  init = 0;
end

if isempty(st) | init
  st = 1;
end

evParms = [];
evParms.SD = iqSet.evParms.SD;
evParms.state.semaphoreOpenedState = 0;
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
evParms.SD.forceInit = 0;
evParms.flag.forceRecord=0;
assignin('base', 'evParms', evParms);
%disp(['nrb = ' num2str(evParms.SD.narrowbandFlag)])

sgp = feval(@spectralDoppler_4p9p2_saveiq,'spectrogramParameters', ...
           evParms.ev.framePeriod_s,...
           evParms.SD.dopSweepTimeIndex);

sweepTime = sgp.THist;

%st=1;
en=st+nPulses-1;

%lastPulseNo = size(IQMat, 3);
lastPulseNo = iqSet.blk_end;


nDopPRIs = iqSet.head.nDopPRIs;
framePeriod = evParms.ev.framePeriod_s;
dopPRF = iqSet.head.dopPRF;
dopWFCutoffNorm  = evParms.SD.dopWFCutoffNorm;
dopSweepTimeIndex = evParms.SD.dopSweepTimeIndex;
baseLineShiftNorm = evParms.SD.baseLineShiftNorm;
SDDynRangeDB = evParms.SD.SDDynRangeDB;
NWindow = evParms.SD.NWindow;
SDDespeckle = evParms.SD.despeckle;


rSel = r; %85:87;
cSel = c; %85:87
%keyboard
while en <= lastPulseNo

  % mode can be switched by callback in offlinesdlargecine
  interactiveMode = evalin('base', 'interactiveMode');

  %disp(['interactiveMode = ' num2str(interactiveMode)])

  if interactiveMode
    r = evalin('base', 'r');
    rSel = [r];    
    c = evalin('base', 'c');  
    cSel = c;
    evalin('base', 'evParms.SD.narrowbandFlag = 0;'); % don't use narrowbanding
    if ishandle(iqSet.evParms.SD.hplotTopInd)
      set(iqSet.evParms.SD.hplotTopInd, 'color', 'b')
    end
  else
    evalin('base', 'evParms.SD.narrowbandFlag = 1;'); % find mean of multiple gates
    evalin('base', 'evParms.SD.narrowbandIndFlag = 1;');
    if ishandle(iqSet.evParms.SD.hplotTopInd)
      set(iqSet.evParms.SD.hplotTopInd, 'color', 'g')
    end
    %    disp('narrowband mode')    
  end
    
  inIQ = IQMat(rSel, cSel, st:en);
  %  size(inIQ)
  %  [rSel, cSel]
  % this is called on each nPulses ensemble
  spectralDoppler_4p9p2_saveiq(real(inIQ), imag(inIQ));
  %pause(50e-3)
  
  stopFlag = evalin('base', 'stopFlag');
  
  if stopFlag
    break
  end
  
  %keyboard
  outStructExists = evalin('base', 'exist(''sdOut'', ''var'')');
  if outStructExists
    assignin('base', 'SDop', SDop);
    evalin('base', 'sdOut.SDop = SDop;');
  end

  tStart_s = (st-1)/evParms.ev.dopPRF;
  tEnd_s = (en-1)/evParms.ev.dopPRF;  
  
  %  title(gca(SDop.hFigure.Children.Children(1)),['start time: ' num2str(tStart_s) ' s'])
  
  gca_this = SDop.hFigure.Children.Children(1);

  %  keyboard
  %sweepTime_s = sgp.THist;
  xt = get(gca_this, 'xtick');
  %   keyboard
  %xtl = num2cell(xt+tStart_s);
  %set(gca_this, 'xticklabel', cell2str(xtl))
  set(gca_this, 'xticklabel', split(sprintf('%1.1f|',xt+tStart_s), '|'))

  gca_BCDI = SDop.hFigure.Children.Children(2);
  %bcdi_tAx_s = evalin('base', 'tAx_s');

  %  ind = find(bcdi_tAx_s >= tStart_s, 1);
  d2s = 24*3600;
  seg_time_s = (iqSet.set_start_datenum-iqSet.set_start_datenum(1))*d2s;
  
  ind = find(seg_time_s <= tEnd_s, 1, 'last');
  hTile = evalin('base', 'hTile');

  set(hTile.imBCDI, 'CData', evalin('base', ['imCompNorm(:,:,:,' num2str(ind) ');']));

  [dum, title_pre] = fileparts(iqSet.matInPath);
  %  inFilePre = split(iqSet.inFilePre, '_sg');

  tit = [];
  tit{1} =  titstrfn(['Set: ' title_pre ':' iqSet.inFilesPre{ind}]);
  tit{2} =  {[ ' at start of segment: ' num2str(ind)]};
  title(hTile.tile, [tit{1} tit{2}]);

  restartFlag = evalin('base', 'restartFlag');

  if restartFlag
    en = 0;
    restartFlag = 0;
    assignin('base', 'restartFlag', restartFlag);
  end
      
  st=en+1;
  en=en+nPulses;

%  pausede
end

% make persistent st start at beginning
if en > lastPulseNo
  st = 1;
end
      

  



end

