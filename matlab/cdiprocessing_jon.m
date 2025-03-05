function IQOut = cdiprocessing_jon(IQData);
% Notice: 
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility 
%   for its use.
%
%cdiprocessing: VDAS external process object function
%John Flynn April 2012

persistent frameCount startTime
persistent pavgHist
persistent Npri NstartPri    FrameHistory
persistent khist
persistent WF H hh
persistent dateStr
persistent fid

numFramesSave=100;
HistoryLength = 100;
if isempty(pavgHist),
    [pavgHist,pimHist,adaptMetricHist]=deal(zeros(1,HistoryLength));
end

if isempty(khist),khist = 0;end
khist = mod(khist,HistoryLength)+1;

if isempty(frameCount)
    frameCount = 0;
    startTime = clock;
    frameRate = 0;
    dateStr = [];
    
    assignin('base', 'IQDataSize', [size(IQData) numFramesSave]);
    assignin('base', 'IQDataSample', IQData(1));
    
    evalin('base', ['IQSave = zeros(IQDataSize, ''like'', ' ...
                    'IQDataSample);']);
    evalin('base', ['IQFrameTime = zeros(' num2str(numFramesSave) ',1);']);

else
    if rem(frameCount,200)==0 ,
        RunningTime = etime(clock,startTime);
        frameRate = frameCount/RunningTime;
        disp([mfilename,':framerate (fps): ',num2str(frameRate)])
        %  size(IQData)
    end
end
frameCount = frameCount+1;
if 1
  frameDst = mod(frameCount-1, numFramesSave)+1;
  assignin('base', 'frameDst', frameDst);
  dimIQ = length(size(IQData));
  assignin('base', 'IQData', IQData);
  %evalin('base', 'size(IQSave)');
  %evalin('base', 'frameDst');
  evalin('base', 'IQSave(:,:,:,:,frameDst) = IQData;');
  evalin('base', ['IQFrameTime(frameDst) = ' num2str(etime(clock,startTime)) ...
                  ';']);
  if frameDst==1
    disp('IQ loop storage started');
  end
end

%assignin('base', 'IQSave(:,:,
if 0
  disp([mfilename,': save frame:'])
  vtHome = getenv('VITTAMED_HOME');
  vtData = [vtHome '_data'];
  matOutPath = [vtData '/iqdata/'];
  if isempty(dateStr)
    fmt=30; %datetime format
    dateStr = datestr(now,fmt)
  end
  
  dateFileName = [ matOutPath mfilename  dateStr '_fr' num2str(frameCount) '.mat'];           
  save(dateFileName, 'IQData');
  lslrt(dateFileName);
end

if 0
  disp([mfilename,': save frame:'])
  vtHome = getenv('VITTAMED_HOME');
  vtData = [vtHome '_data'];
  matOutPath = [vtData '/iqdata/'];
  if isempty(dateStr)
    fmt=30; %datetime format
    dateStr = datestr(now,fmt)
  end
  
  dateFileName = [ matOutPath mfilename  dateStr '_fr' num2str(frameCount) '.mat'];           
  save(dateFileName, 'IQData');
  lslrt(dateFileName);
end
             
%wall filtering, AC , detection, and normalize
%this method replaces IQOut at NstartPri with AC result.
persistent NfiltBank

Process = evalin('base','Process');
procIndCFI = evalin('base',  'procIndCFI');
P = evalin('base','P');
PIndCFI = evalin('base','PIndCFI');

pc = {Process.classname};

%pwrThres = evalin('base','pwrThres');
pwrThres = P(PIndCFI).pwrThres;
DopState = P(PIndCFI).DopState;
imgPID =  procIndCFI; %must match script, doppler Imag process % strmatch('Image',pc,'exact') ;

DopState = evalin('base','DopState');
    imgMeth = 'imageDisplay';
 
switch (DopState)
    case     'computeCFIPowerEst';
        PowerEstMode = 1;
    case      'computeCFIFreqEst';
        PowerEstMode = 0;
    otherwise
        error('bad switch')
end

[Nr, Nc, Nf, Np  ]=size(IQData);

%compute and store wall filter coefficients
if isempty( WF ) ,
    NstartPri = 3; %throw away beginning pulses in ensemble
    Npri = Np - NstartPri + 1;
    ord = ceil(Npri/8);
    h  = polyfilter(Npri,ord);
    
    H = single(h');%note: RIGHT HAND MULTIPLY
    [WF ] = single(zeros(Np)) ;
    WF(NstartPri:end,NstartPri:end)= H  ;
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing

%assume only one frame to process with WF
IQD = reshape(   IQData   , Nr*Nc*Nf , Np );

%wallfiltering
IQDwf = IQD*WF;

%lag-1 autocorrelation
ac1 = mean(IQDwf(:,2:end).*conj(IQDwf(:,1:end-1)),2);
ac1=ac1.';

% detect and normalize:
threshParams = struct('threshold', pwrThres ,'PowerEstMode',PowerEstMode,'Npri',Npri);
[img, pow, threshData] = vsCDIAdaptiveThreshold(  ac1 , threshParams );

IQOut = reshape(  double(img)   , Nr, Nc );

end %main %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

