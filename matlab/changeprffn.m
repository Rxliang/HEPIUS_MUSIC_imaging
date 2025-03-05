function Control = changeprffn(UIValue, hObject, deferControlUpdate)

if nargin < 3
  deferControlUpdate=0;
end

% UI value is the Doppler PRF.
% For interleaved B/SD, ttna is 1/(2*UIValue/1e-6)
% For SD only, it is 1/(UIValue/1e-6)

evParms = evalin('base', 'evParms');
evParmsOrig = evParms;

Resource = evalin('base', 'Resource');
ResourceOrig = Resource;

Receive = evalin('base', 'Receive');
ReceiveOrig = Receive;

evParms.state.SDOnly=0;
SDOnly= evParms.state.SDOnly;

dopPRFNom = UIValue;
framePeriod_s = evParms.ev.framePeriod_s;

if evParms.state.SDOnlyInPlace
  ttna_microsec = round(1/dopPRFNom/1e-6);
  nPRIs = round(framePeriod_s/ttna_microsec*1e6);
else
  ttna_microsec = round(1/dopPRFNom/2e-6);
  nPRIs = round(1/2*framePeriod_s/ttna_microsec*1e6)*2; % make it even
end

if nPRIs > evParms.ev.maxPRIs
  disp('nPRIs > max. setting to max');
  nPRIs = evParms.ev.maxPRIs
  ttna_microsec = round(1e6*framePeriod_s/nPRIs)
end

evParms.ev.ttna_microsec = ttna_microsec;
evParms.ev.nPRIs = nPRIs;

if evParms.state.SDOnlyInPlace % remove interleaving for SD only mode
  evParms.ev.dopPRF=1/(evParms.ev.ttna_microsec*1e-6);
  nPRIs = round(framePeriod_s/ttna_microsec*1e6);
  evParms.ev.priSkip = 1;
  evParms.state.SDOnly = SDOnly;
  evParms.ev.nDopPRIs = nPRIs;
else
  evParms.ev.dopPRF=1/(evParms.ev.ttna_microsec*2e-6);
  evParms.ev.nDopPRIs = evParms.ev.nPRIs/2;  
  evParms.ev.priSkip = 2;
end

SeqControlOrig = evalin('base', 'SeqControl');
ReceiveOrigIn=[];

seqContainer = makenicpseqfn(evParms, Resource, SeqControlOrig, ...
                             ReceiveOrigIn);

evParms = seqContainer.evParms;
recon = seqContainer.recon;
ri = seqContainer.ri;

Receive = seqContainer.Receive;
Recon = seqContainer.Recon;
ReconInfo = seqContainer.ReconInfo;
Resource = seqContainer.Resource;
Process = seqContainer.Process;
Event = seqContainer.Event;
SeqControl = seqContainer.SeqControl;

% if we increase PRF, new Events will reference Receives that do
% not exist. So when reducing PRF, keep these. don't help. fails
% VSX check

%if length(Receive) < length(ReceiveOrig)
%  Receive(length(Receive)+1:length(ReceiveOrig)) = ReceiveOrig(length(Receive)+1: ...
%                                                 end);
%end

ReconOrig= evalin('base', 'Recon');
ReconInfoOrig = evalin('base', 'ReconInfo');

if 1 % this is necessary. If it is removed, will get crash since
     % numchannels will be missing
      
if isfield(ReconOrig, 'numchannels') % & (length(Recon) < length(ReconOrig)  )
  for q = 1:length(Recon)
    Recon(q).numchannels = ReconOrig(q).numchannels;
  end
end

if isfield(ReconInfoOrig, 'Aperture') % & (length(ReconInfo) < length(ReconInfoOrig))
  for q = 1:length(ReconInfo)
    ReconInfo(q).Aperture = ReconInfoOrig(1).Aperture;
  end
end

end

l1=length(ReconInfo); l0 = length(ReconInfoOrig);
if l1 < l0
  ReconInfo(l1+1:l0)=ReconInfoOrig(l1+1:l0);
end

if 0
addInd = 2*nPRIs+1:3*nPRIs;
ReconInfo(addInd) = ReconInfo(1:nPRIs);
Recon(1).RINums = addInd;
%Recon(2).RINums = 251:375;
end

% this is used if we reduce number of events below those referenced
% by SeqControl
%SeqControl = resetseqctrljumpfn(SeqControl);
 
EventOrig = evalin('base', 'Event');
ProcessOrig = evalin('base', 'Process');

lenEventOrig = length(EventOrig);

lenEvent= length(Event);

if 0
  % this idea of padding does not help avoid: Error using runAcq
  % getStruct: Event.seqControl with multiple values can't have any
  % values set to 0.

  % make sure we never reduce number of events below max
  persistent maxEvents
  if isempty(maxEvents)
    maxEvents = lenEventOrig;
  end

  maxEvents = max(lenEvent, maxEvents);
  
  if lenEvent < maxEvents
    % pad with useless events
%    for q = lenEvent+1:evParms.ev.numEventMax
    for q = lenEvent+1:maxEvents
      Event(q) = eventnoopfn(evParms.seq);
    end
  end
end

%Event(lenEvent+1:lenEventOrig) = EventOrig(lenEvent+1:lenEventOrig);

if 1
  disp('assignin');
%  pausede
% reinitialize sdfn
%evParms.SD.forceInit = [1 1];

assignin('base', 'evParms', evParms);
assignin('base', 'recon', recon);
assignin('base','ri', ri);
assignin('base','Receive', Receive);
assignin('base', 'Recon', Recon);
assignin('base', 'ReconInfo', ReconInfo);
% PRF change will change CFI processes' PRF parameter
assignin('base', 'Process', Process);
%assignin('base', 'Resource', Resource);
assignin('base', 'Event', Event);
assignin('base', 'SeqControl', SeqControl);
end


if 0
diffstructarrfn(evParms, evParmsOrig);
diffstructarrfn(Receive, ReceiveOrig);
diffstructarrfn(Recon, ReconOrig)
diffstructarrfn(Event, EventOrig)
diffstructarrfn(SeqControl, SeqControlOrig);
diffstructarrfn(ReconInfo, ReconInfoOrig);
diffstructarrfn(Process, ProcessOrig)
diffstructarrfn(Resource, ResourceOrig)
end

if 1
if evalin('base', 'exist(''VDASupdates'', ''var'')')
  disp('vsupdate run');
  evalin('base','VsUpdate(''Receive'')');
%  evalin('base','VsUpdate(''SeqControl'')');  
end
end

Control = evalin('base','Control');

ctrlEmpty = pendingcontrolprocfn(evParms.lock.ctrlLockSleepTime_s);

if ctrlEmpty
  c=1;
else
  c=length(Control)+1;
end

Control(c).Command = 'set&Run';
Control(c).Parameters = {'Parameters',1,'startEvent', 1}; 
evalin('base',['Resource.Parameters.startEvent = ' ...
               num2str(1) ';']);  

if 1
c=c+1;
Control(c).Command = 'set&Run';
Control(c).Parameters = {'Process',evParms.proc.procIndB, ...
                         'pgain', evParms.ev.BImGain};

c=c+1;
if evParms.state.SDOnlyInPlace
  % SDLarge keeps using the 2-spaced PRIs, so it has half the
  % SD freq
  procSDLargePRF = evParms.ev.dopPRF/2;
else
  procSDLargePRF = evParms.ev.dopPRF/evParms.largeSDParms.priSkip;
end
Control(c).Command = 'set&Run';
Control(c).Parameters = {'Process',evParms.proc.procIndSDLarge, ...
                          'prf',  procSDLargePRF};

end

c=c+1;
Control(c).Command = 'update&Run';
Control(c).Parameters = {'Recon', 'Event','SeqControl','Receive'};

if ~deferControlUpdate
  assignin('base', 'Control', Control);
end

if nargin > 1 & ~isempty(hObject) & ishandle(hObject)
  set(hObject,'Value',  evParms.ev.dopPRF);
end

% Update the WF control:
hwf = findobj('tag','UserB3Slider');
dopWFCutoffNorm = evParms.SD.dopWFCutoffNorm;
set(hwf,'Value',evParms.ev.dopPRF*dopWFCutoffNorm);
hwf = findobj('tag','UserB3Edit');
set(hwf,'String',num2str(evParms.ev.dopPRF*dopWFCutoffNorm,'%4.0f Hz'));
% update PRF slider display
hwf = findobj('tag','UserB4Edit');
set(hwf,'String',num2str(evParms.ev.dopPRF,'%4.0f Hz'));








