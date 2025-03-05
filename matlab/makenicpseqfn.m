function seqContainer = makenicpseqfn(evParms, Resource, SeqControl, ...
                                      Receive, seq)

if nargin < 3
  SeqControl = [];
  seq = [];
else
  if isfield(evParms, 'seq')
    seq = evParms.seq;
  else
    seq = [];
  end
  
end

defaultLoopTstArg = 1; % placeholder event to jump to in case
                       % SDOnly has fewer events
defaultJumpArg = 1; % placeholder event to jump to in case

rcvParms = evParms.rcvParms;

% both needed so we can keep number of SeqControl equal prospectively

if evParms.state.SDOnly
  evParms.flag.doCFINow=0;
end

if evParms.state.SDOnlyInPlace
  evParms.ev.nDopPRIs=evParms.ev.nPRIs;
  evParms.ev.dopPRF=evParms.ev.acqPRF;
  evParms.ev.priSkip=2; % each loop adds two Events
end

ReceiveOrig=Receive;
[all.Receive, all.rcvOut] = sdcallreceivefn(rcvParms, evParms, ...
                                         evParms.P, ReceiveOrig);

if ~evParms.state.SDOnly
  Receive = all.Receive;
else
  Receive = all.Receive;
end

rcvOut = all.rcvOut;
rcvOut.PIndCFI = rcvParms.PIndCFI;
evParms.rcvOut = rcvOut;

Recon = []; ri = []; recon = [];

if ~evParms.state.SDOnly
  % this makes all the recons: B, SD, CFI
  [bsd.Recon, bsd.ReconInfo, bsd.recon, bsd.ri] = sdcreconfn(Recon, recon, ...
                                                  ri, evParms, ...
                                                  rcvParms);
else
  uscmsg(evParms.state.execState);
  [bsd.Recon, bsd.ReconInfo, bsd.recon, bsd.ri] = nicpsdonlyreconfn(recon, ...
                                                    ri, evParms, rcvParms);
end


% this has all we need for B/SD/CFI
Recon=bsd.Recon;
ReconInfo=bsd.ReconInfo;
recon=bsd.recon;
ri = bsd.ri;

if ~evParms.state.SDOnly
  [Process, proc, Resource, evParms] = sdcprocfn(evParms, Resource);
else
  uscmsg(evParms.state.execState);
  [Process, proc, Resource, evParms] = nicpsdonlyprocfn(evParms, Resource);
end


if evParms.flag.runROIPlotForSDOnly
  sdo.evIndStartAcq = 3;
else
  sdo.evIndStartAcq = 2;
end

bsd.evIndStartB = 3; % will check this is correct later
bsd.evIndStartAcq = 6; % start of acq loop
cfi.evIndStartB = 2; % will check this is correct later
cfi.evIndStartAcq = 5; % start of acq loop

all.evIndStartB = 3; % will check this is correct later
all.evIndStartAcq = 6; % start of acq loop

nsc=1;

evIndStartAcq=all.evIndStartAcq;
evIndStartB = all.evIndStartB;

SeqControl = [];
%if ~evParms.state.SDOnly
  keepSeqCtrlNonSDOnly=1;
%else
 % keepSeqCtrlNonSDOnly=0;
%end

if isempty(SeqControl)
  if keepSeqCtrlNonSDOnly
  % this is never used explicitly, but VSX wants a SeqControl(1)
  SeqControl(nsc).command = 'jump';
  SeqControl(nsc).argument = evIndStartB; % start of 2D acq
  seq.seqControlIndJumpToB = nsc;
  nsc=nsc+1;
  end
  
  SeqControl(nsc).command = 'timeToNextAcq';
  SeqControl(nsc).argument = evParms.ev.ttna_microsec ; % time in usec
  seq.seqControlIndTTNASD = nsc;
  % loopCnt/loopTst are used to implement a conditional jump
  nsc=nsc+1;
  
  if keepSeqCtrlNonSDOnly
  % set loopCnt to zero
  
  SeqControl(nsc).command = 'loopCnt';
  SeqControl(nsc).argument = 0;
  SeqControl(nsc).condition = 'counter1';
  seq.seqControlIndSetLoopCnt0 = nsc;

  % set loopCnt to large value
  nsc=nsc+1;
  SeqControl(nsc).command = 'loopCnt';
  % allow 1 jump per frame
  SeqControl(nsc).argument = evParms.ev.nfrms; % can be up to 65536
  SeqControl(nsc).condition = 'counter1';
  seq.seqControlIndSetLoopCnt1 = nsc;

  nsc=nsc+1;
  SeqControl(nsc).command = 'jump';
  SeqControl(nsc).argument = evIndStartAcq;  % jump to start of acquisition
  seq.seqControlIndJumpStartEv = nsc;
  nsc=nsc+1;

  % -- Change to TPC Profile SD
  SeqControl(nsc).command = 'setTPCProfile';
  SeqControl(nsc).condition = 'next';
  SeqControl(nsc).argument = evParms.ind.TPCIndSD;
  seq.seqControlIndSetTPCProfileSD= nsc;
  nsc = nsc + 1;

  SeqControl(nsc).command = 'setRcvProfile';        
  SeqControl(nsc).argument = evParms.ind.RcvProfileIndSD;  
  seq.seqControlIndSetRcvProfileSD = nsc;
  nsc = nsc + 1;

  % end of frame time for multiangle plane wave
  SeqControl(nsc).command = 'timeToNextAcq';
  SeqControl(nsc).argument = evParms.ev.transferTTNA_us;
  seq.seqControlIndTransferTTNA = nsc;
  nsc = nsc + 1;

  % end of frame time for multiangle plane wave
  SeqControl(nsc).command = 'timeToNextAcq';
  SeqControl(nsc).argument = evParms.ev.transferTTNA_postSD_us;
  seq.seqControlIndTransferTTNA_postSD = nsc;
  nsc = nsc + 1;
  
  % these are needed for evParms.flag.doCFI, but since we cannot
  % change number of SeqControl elements, define whether needed or
  % not
  % - set receive profile for 2D
  SeqControl(nsc).command = 'setRcvProfile';        
  SeqControl(nsc).argument = evParms.ind.RcvProfileIndB;  
  seq.seqControlIndSetRcvProfile2D = nsc;
  nsc = nsc + 1;

  % - time between wide beams 
  SeqControl(nsc).command = 'timeToNextAcq';
  SeqControl(nsc).argument = 200;
  nsc = nsc + 1;
  
  % -- Time between 2D+SD acquisition and Doppler ensemble. 
  % Set to allow time for profile change.
  SeqControl(nsc).command = 'timeToNextAcq';
  SeqControl(nsc).argument = 2000; % time in usec
  seq.seqControlIndSetLongTTNACFI = nsc;
  seq.seqControlIndSetLongTTNAChangeProfile = nsc;
  nsc = nsc + 1;

  % -- Change to TPC Profile CFI (Doppler CFI)
  SeqControl(nsc).command = 'setTPCProfile';
  SeqControl(nsc).condition = 'next';
  SeqControl(nsc).argument = evParms.ind.TPCIndCFI;
  seq.seqControlIndSetTPCProfileCFI = nsc;
  nsc = nsc + 1;
  
  % -- Change to TPC Profile CFI (Doppler CFI)
  SeqControl(nsc).command = 'setTPCProfile';
  SeqControl(nsc).condition = 'immediate';
  SeqControl(nsc).argument = evParms.ind.TPCIndCFI;
  seq.seqControlIndSetTPCProfileCFIImmed = nsc;
  nsc = nsc + 1;
  
  % -- Change to TPC Profile SD immed
  SeqControl(nsc).command = 'setTPCProfile';
  SeqControl(nsc).condition = 'immediate';
  SeqControl(nsc).argument = evParms.ind.TPCIndSD;
  seq.seqControlIndSetTPCProfileSDImmed = nsc;
  nsc = nsc + 1;
  end
  
  % - set Receive profile for Doppler
  SeqControl(nsc).command = 'setRcvProfile';
  SeqControl(nsc).argument =  evParms.ind.RcvProfileIndCFI;
  seq.seqControlIndSetRcvProfileCFI = nsc;
  nsc = nsc + 1;

  if keepSeqCtrlNonSDOnly
  % - time between Doppler ensemble elements (P(PIndCFI).dopPRF)
  SeqControl(nsc).command = 'timeToNextAcq';
  seq.seqControlIndSetShortTTNACFI = nsc;
  % argument is set outside this SeqControl initialization conditional block
  nsc = nsc + 1;
  
  % not needed for no CFI case, but define anyway for that user
  % switch case
   % - set Receive profile for Doppler, even though not doing CFI
  SeqControl(nsc).command = 'setRcvProfile';
  SeqControl(nsc).argument =  evParms.ind.RcvProfileIndCFI;
  seq.seqControlIndSetRcvProfileCFI = nsc;
  nsc = nsc + 1;

  % -- Change to Profile 1 (2D)
  SeqControl(nsc).command = 'setTPCProfile';
  SeqControl(nsc).condition = 'next';
  SeqControl(nsc).argument = evParms.ind.TPCIndB;
  seq.seqControlIndSetTPCProfileB = nsc;
  nsc = nsc + 1;

  % - time between frames  
  SeqControl(nsc).command = 'timeToNextAcq';
  SeqControl(nsc).argument = 2000;
  nsc = nsc + 1;
  end
  
  % - return to Matlab
  SeqControl(nsc).command = 'returnToMatlab';
  SeqControl(nsc).argument = 0;  
  seq.seqControlIndReturnToMatlab = nsc;

  if evParms.flag.sendTriggers
    nsc=nsc+1; % for debugging, produce a hardware trigger out
               % set trigger out
    SeqControl(nsc).command = 'triggerOut';        
    SeqControl(nsc).argument = 0;  % no delay 
    seq.seqControlIndTriggerOut = nsc;
  end
  nsc=nsc+1;
  
  if keepSeqCtrlNonSDOnly
  SeqControl(nsc).command = 'sync';
  SeqControl(nsc).argument = 250000; % us
  seq.seqControlIndSync=nsc;
  nsc=nsc+1;

  noopDelay_s = evParms.ev.preFrameHWDelay_s;
  % in units of 200ns
  noopDelayArg = round(noopDelay_s/200e-9);
  SeqControl(nsc).command = 'noop';
  SeqControl(nsc).argument = noopDelayArg; 
  seq.seqControlIndNoop=nsc;
  nsc=nsc+1;
  end
  
  % preallocate transfers so we don't mess with these between
  % sequence swicthes
  seq.seqControlIndLoopTstSkipCFI = [];
  seq.seqControlIndTransferToHost = [];
  seq.seqControlIndJumpOverCFIAcq = [];
  seq.seqControlIndJumpOverSDRecon = [];
  seq.seqControlIndJumpOverCFIRecon = [];
  
  for q = 1:evParms.ev.nfrms
    if keepSeqCtrlNonSDOnly
    % loopTst
    SeqControl(nsc).command = 'loopTst';
    SeqControl(nsc).condition = 'counter1';    
    SeqControl(nsc).argument = defaultLoopTstArg;
    seq.seqControlIndLoopTstSkipCFI(q) = nsc; 
    nsc=nsc+1;
    end
  
    SeqControl(nsc).command = 'transferToHost';
    seq.seqControlIndTransferToHost(q) = nsc;
    nsc=nsc+1;
    
    if keepSeqCtrlNonSDOnly
    % jump over CFI acq
    SeqControl(nsc).command = 'jump'; 
    SeqControl(nsc).argument = defaultJumpArg;
    seq.seqControlIndJumpOverCFIAcq(q) = nsc; 
    nsc=nsc+1;

    % loopTst jump over SD recon
    SeqControl(nsc).command = 'loopTst';
    SeqControl(nsc).condition = 'counter1';    
    SeqControl(nsc).argument = defaultLoopTstArg;
    seq.seqControlIndJumpOverSDRecon(q) = nsc;
    nsc=nsc+1;

    % jump over CFI recon
    SeqControl(nsc).command = 'jump';
    SeqControl(nsc).argument = defaultJumpArg;
    seq.seqControlIndJumpOverCFIRecon(q) = nsc;  
    nsc=nsc+1;
    end
    
  end
    
  % Jump back to first event
  SeqControl(nsc).command = 'jump';
  SeqControl(nsc).argument = defaultJumpArg;
  seq.seqControlIndJumpToFirstEvent = nsc;
  nsc=nsc+1;
  
end % initialization of SeqControl block

% SD PRF adjustment
SeqControl(seq.seqControlIndTTNASD).argument = evParms.ev.ttna_microsec; 
    
% conditional adjustment of SeqControl argument

if keepSeqCtrlNonSDOnly
if evParms.flag.doCFISeq
  SeqControl(seq.seqControlIndSetShortTTNACFI).argument = ...
      round(1/(P(evParms.ind.PIndCFI).dopPRF*1e-06));
else
  SeqControl(seq.seqControlIndSetShortTTNACFI).argument = ...
      round(1/(evParms.P( evParms.ind.PIndCFI).dopPRF*1e-06))...
      / evParms.ev.priSkip;
end
end

% this is problematic, since I don't think we can change
% seqControl on the fly. may need to do something else here ***
 
if ~evParms.flag.doCFISep
  %---------------check Doppler PRF--------------------
  tof = ceil(2*Receive(rcvCFIStartInd(1)).endDepth/demodFreq);
  if SeqControl(seq.seqControlIndSetShortTTNACFI).argument < tof
    SeqControl(seq.seqControlIndSetShortTTNACFI).argument = tof;
    P(evParms.ind.PIndCFI).dopPRF = round(1/(tof*1e-06));
    SeqControl(seq.seqControlIndSetShortTTNACFI).argument = ...
        round(1/(P(evParms.ind.PIndCFI).dopPRF*1e-06));
    fprintf(['"timeToNextAcq" is adjusted to ' num2str(tof) '\n']);
    fprintf(['"P.dopPRF" is adjusted to ' num2str(P(evParms.ind.PIndCFI).dopPRF) '\n']);
    for k = 1:2:length(Process(proc.procIndCFI).Parameters)
      if strcmp(Process(proc.procIndCFI).Parameters{k},'prf'), ...
            Process(proc.procIndCFI).Parameters{k+1} = P(evParms.ind.PIndCFI).dopPRF;
      end       
    end
  end
end


evParms.ev.evIndStartB = evIndStartB;
evParms.ev.evIndStartAcq = evIndStartAcq;
evParms.proc = proc;
evParms.seq = seq;
evParms.recon = recon;
evParms.SeqControlPreEv = SeqControl;

cfi.evParms=evParms;
cfi.evParms.ev.evIndStartB = cfi.evIndStartB;
cfi.evParms.ev.evIndStartAcq = cfi.evIndStartAcq;
%cfi.evParms.ev.TXIndCFI = TXIndCFI;
%cfi.evParms.ev.rcvCFIStartInd = rcvCFIStartInd;

bsd.evParms=evParms;
bsd.evParms.ev.evIndStartB = bsd.evIndStartB;
bsd.evParms.ev.evIndStartAcq = bsd.evIndStartAcq;
bsd.evParms.ev.TXOffset =  evParms.ev.TXIndCFI-1;

all.evParms=evParms;
all.evParms.ev.evIndStartB = all.evIndStartB;
all.evParms.ev.evIndStartAcq = all.evIndStartAcq;
all.evParms.ev.TXOffset = evParms.ev.TXIndCFI-1;

sdo.evParms=evParms;
sdo.evParms.ev.evIndStartAcq = sdo.evIndStartAcq;

evParms = all.evParms;

%if ~evParms.state.SDOnly
%  cfi.evParms.rcvOut=cfi.rcvOut;
%end

%bsd.evParms.rcvOut=bsd.rcvOut;
%bsd.evParms.recon=bsd.recon;

% save the non-sequence dependent SeqControl for later use as
% starting point

if ~evParms.state.SDOnly
  
  [all.Event, all.SeqControl, all.evParmsOut] = sdcalleventfn(all.evParms, ...
                                                    SeqControl);

  Event = all.Event;
  SeqControl =  all.SeqControl;
  evParms = all.evParmsOut;

else
  uscmsg(evParms.state.execState);
  [sdo.Event, sdo.SeqControl, sdo.evParmsOut] = ...
      sdceventsdonlyfn(sdo.evParms, SeqControl);
  Event = sdo.Event;
  SeqControl =  sdo.SeqControl;
  evParms = sdo.evParmsOut;
end

seqContainer.evParms = evParms;
seqContainer.recon = recon;
seqContainer.ri = ri;
% VSX vars
seqContainer.Receive = Receive;
seqContainer.Recon = Recon;
seqContainer.ReconInfo = ReconInfo;
seqContainer.Resource = Resource;
seqContainer.Process = Process;
seqContainer.Event = Event;
seqContainer.SeqControl = SeqControl;


return

if 0
if ~evParms.state.SDOnly
  if evParms.flag.doCFISepAtStart
    Event = cfi.Event;
    SeqControl = cfi.SeqControl;
    evParmsOut = cfi.evParmsOut;
    
    if length(bsd.SeqControl) >  length(cfi.SeqControl)
      SeqControl(length(bsd.SeqControl)).command  = [];
    end
%    [Event, SeqControl, evParmsOut] = sdccfieventfn(evParms, SeqControl);
  else
    Event = bsd.Event;
    SeqControl = bsd.SeqControl;
    evParmsOut = bsd.evParmsOut;
%    [Event, SeqControl, evParmsOut] = sdceventfn(evParms, ...
%                                                 SeqControl);
    if length(cfi.SeqControl) >  length(bsd.SeqControl)
      SeqControl(length(bsd.SeqControl)).command  = [];
    end
    
      
  end
  evParms=evParmsOut; 
  evParms.seq.lenSeqControl = length(SeqControl);
else
  dopPRF=evParms.ev.acqPRF;
%  ttna_microsec = psched.PRIAcq_uSec;
  evParms.ev.nDopPRIs = evParms.ev.nPRIs;
  evParms.ev.priSkip = 1;
  evParms.ev.evIndStartAcq = 2;
  [Event, SeqControl, evParmsOut] = sdceventsdonlyfn(evParms);
%  [Event, SeqControl, evParmsOut] = sdceventbonlyfn(evParms);
end
end

end

