function [IQOut, rawOut] = cdiprocessing_nxtv1p0p4(IData, QData)
%
% cdiprocessing: VDAS external process object function
%
% The wall filter coefficients computed by this function correspond
% to polynomial regression filters having stop-band order commensurate
% with ensemble length. This differs from the runAcq "regression" wall
% filter type, which has a fixed stop-band order.
%
% Notice:
%   This file is provided by Verasonics to end users as a programming
%   example for the Verasonics Vantage NXT Research Ultrasound System.
%   Verasonics makes no claims as to the functionality or intended
%   application of this program and the user assumes all responsibility
%   for its use.
%
% Copyright © 2013-2023 Verasonics, Inc.

persistent frameCount startTime
persistent pavgHist
persistent Npri NstartPri    FrameHistory
persistent khist
persistent WF H hh

IQData = complex(IData, QData); % create complex array used by this function

HistoryLength = 100;
if isempty(pavgHist),
    [pavgHist,pimHist,adaptMetricHist]=deal(zeros(1,HistoryLength));
end

if isempty(khist),khist = 0;end
khist = mod(khist,HistoryLength)+1;

if isempty(frameCount)
    frameCount = 0;
    startTime =clock;
    frameRate = 0;
else
    if rem(frameCount,200)==0 ,
        RunningTime = etime(clock,startTime);
        frameRate = frameCount/RunningTime;
        disp([mfilename,':framerate (fps): ',num2str(frameRate)])
        %  size(IQData)
    end
end
frameCount = frameCount+1;

%wall filtering, AC , detection, and normalize
%this method replaces IQOut at NstartPri with AC result.
persistent NfiltBank

%Process = evalin('base','Process'); %JSM
%pc = {Process.classname};

pwrThres = evalin('base','pwrThres');
% JSM imgPID =  3; %must match script, doppler Imag process % strmatch('Image',pc,'exact') ;

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

% Processing

%assume only one frame to process with WF
IQD = reshape(   IQData   , Nr*Nc*Nf , Np );

%keyboard
%wallfiltering
IQDwf = IQD*WF;

%lag-1 autocorrelation
ac1 = mean(IQDwf(:,2:end).*conj(IQDwf(:,1:end-1)),2);
ac1=ac1.';

% detect and normalize:
threshParams = struct('threshold', pwrThres ,'PowerEstMode',PowerEstMode,'Npri',Npri);
[img, pow, threshData] = vsCDIAdaptiveThreshold(ac1, threshParams );

IQOut = reshape(  double(img)   , Nr, Nc );
rawOut = reshape(  double(ac1)   , Nr, Nc );

end

