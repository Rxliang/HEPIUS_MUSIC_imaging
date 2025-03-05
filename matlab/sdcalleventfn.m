function [Event, SeqControl, evParmsOut] = sdcalleventfn(evParms, SeqControl)

zeroB = 0;
zeroBRecon=0;
evParmsOut=evParms;

if zeroB
  evParmsOut.ev.dopPRF=evParms.ev.dopPRF*2;
  evParmsOut.ev.framePeriod_s =  evParms.ev.framePeriod_s/2;
end

doCFI = evParms.flag.doCFI;
doCFISeq = evParms.flag.doCFISeq;
multiAngleB = evParms.flag.multiAngleB;
proc = evParms.proc;
seq = evParms.seq;
evIndStartB = evParms.ev.evIndStartB;
evIndStartAcq = evParms.ev.evIndStartAcq;
preFrameHWDelay_s =  evParms.ev.preFrameHWDelay_s;
sendTriggers = evParms.flag.sendTriggers;
syncFrameStart = evParms.flag.syncFrameStart;
singleRcvProfile = evParms.flag.singleRcvProfile;
rcvFrameStartInd = evParms.rcvOut.rcvFrameStartInd;
trig = evParms.trig;

priSkip = evParms.ev.priSkip;
TXIndBVec = evParms.ev.TXIndBVec;
TXIndB = evParms.ev.TXIndB;
TXIndSD = evParms.ev.TXIndSD;
lastFrameTTNA =  evParms.ev.lastFrameTTNA;
P = evParms.P;
nfrms =  evParms.ev.nfrms;
singleTpcProfile = evParms.flag.singleTpcProfile;
recon = evParms.recon;

nsc = length(SeqControl)+1;

TXIndCFI = evParms.ev.TXIndCFI;
TXOffset = evParms.ev.TXOffset;
PIndCFI = evParms.ind.PIndCFI;
rcvCFIStartInd = evParms.rcvOut.rcvCFIStartInd;

m = P(PIndCFI).dopNumRays;

n = 1; % n is count for next Event

% always run so we can set power slides and move main window
Event(n).info = 'ROI plot';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process= 0;     
Event(n).seqControl = 0; 
Event(n).process = proc.procIndCFIROIPlot; % outline doppler
                                           % region used in
                                           % doCFISeq 

n = n+1;

% set recv profile
Event(n).info = 'Recv profile';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0; % outline doppler region
Event(n).seqControl =  [seq.seqControlIndSetRcvProfileSD ...
                        seq.seqControlIndSetTPCProfileSDImmed ];

n = n+1;

% B-mode sequence start event:
if evIndStartB ~= n
  error(['Assumption that event ' num2str(evIndStartB) ' is start of B-mode is wrong!']);
end

Event(n).info = 'Start here for B/SD processing.';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % no process

if ~evParms.flag.startInBCFI
  Event(n).seqControl = seq.seqControlIndSetLoopCnt0; % set loopCnt
                                                      % to 0
else
  % start in BCFI
  Event(n).seqControl = seq.seqControlIndSetLoopCnt1; % set loopCnt
end

n = n+1;

Event(n).info = 'Jump to start of acquisitions';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % no process
Event(n).seqControl = seq.seqControlIndJumpStartEv; % jump

n = n+1;
Event(n).info = 'Start here for B/CFI processing.';
Event(n).tx = 0;         % no transmit
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no reconstruction
Event(n).process = 0;    % no process
Event(n).seqControl = seq.seqControlIndSetLoopCnt1; % set loopCnt to 1
evParmsOut.ev.evIndStartBCFI=n;

n = n+1;

evIndStartSD = n;
nStartAcq = n;
if evIndStartAcq ~= n
  error(['Assumption that event ' num2str(evIndStartAcq) ...
         ' is start of acq loop is wrong! Current value is ' num2str(n)]);
end

evIndFrameStart = [];
for i = 1:nfrms  % Acquire frames
    evIndFrameStart(i) = n;
    
    Event(n).info = 'test loopCnt to see if we should skip to CFI';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no process
    %Event(n).seqControl = nsc; % test loopCnt
    Event(n).seqControl = seq.seqControlIndLoopTstSkipCFI(i);
    %SeqControl(nsc).command = 'loopTst';
    %SeqControl(nsc).condition = 'counter1';    
    %seqControlIndBCFIAcq = nsc; % need to set argument below
    %nsc=nsc+1;
    
    n=n+1;
        
    if preFrameHWDelay_s
      Event(n).info = 'Delay to let SW catch up';
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = 0;      % no reconstruction
      Event(n).process = 0;    % no processing
      Event(n).seqControl = seq.seqControlIndNoop;
      n=n+1;
    end

    if syncFrameStart | 1
        disp('SYNC FRAME START');
        Event(n).info = 'Sync';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction
        Event(n).process = 0;    % no processing
        Event(n).seqControl = [Event(n).seqControl ...
                            seq.seqControlIndSync];
        n=n+1;        
    end
    
    if sendTriggers & 0
      % this appears not to work
      Event(n).info = 'Send trigger out';
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = 0;      % no reconstruction
      Event(n).process = 0;    % no processing
      Event(n).seqControl = [seq.seqControlIndTriggerOut];

      
      n=n+1;
    end

    k = rcvFrameStartInd(i)-1;
    txCnt = 0;    
    thisFrameRcvB=[];
    
    if ~evParms.flag.separateBAcq
      priSkip=1;
    end
    
    for j = 1:priSkip:evParms.ev.nPRIs % acquisition repeats every 2
                                       % acquisitions
      if evParms.flag.separateBAcq
        txCnt=txCnt+1;
        Event(n).info = '2D acquisition.';
        if multiAngleB
          Event(n).tx = TXIndBVec(txCnt);  % use TX structure 1.
        else
          Event(n).tx = TXIndB;  % use TX structure 1.
        end
        
        Event(n).rcv = k+j;
        thisFrameRcvB(txCnt) = Event(n).rcv;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = seq.seqControlIndTTNASD; % set ttna
        if ismember(i, trig.frameNos) & (...
           (j==1 & sendTriggers & trig.onFirstPulseOfTypeOnly & trig.onB) | ...
           (sendTriggers & trig.onB & ~trig.onFirstPulseOfTypeOnly) | ...
           (sendTriggers & trig.onSD & ~trig.onFirstPulseOfTypeOnly ...
            & evParms.state.SDOnlyInPlace)...
                                        )
          Event(n).seqControl = [Event(n).seqControl ...
                                 seq.seqControlIndTriggerOut];
        end

        if zeroB
          Event(n).tx=0;
          Event(n).rcv=0;
          Event(n).seqControl=0;
        end
        
        if evParms.flag.BLineWT
          uscmsg(evParms.state.execState);
          rayIndBWT=0;
          if ismember(j, evParms.BLineWT.TXIndBWT)
            Event(n).tx =  evParms.ev.TXIndBWTStart+rayIndBWT;
            Event(n).info = 'WT acquisition';
          end
        end
        
        if evParms.state.SDOnlyInPlace
           Event(n).info = 'Doppler acquisition';
           if j==1 
%             Event(n).tx = TXIndB; % experiment that did not work    
             Event(n).tx = TXIndSD;        
           else
             Event(n).tx = TXIndSD;        
           end
           
           Event(n).rcv = k+j;
           Event(n).recon = 0;      % no reconstruction.
           Event(n).process = 0;    % no processing
           
           % inherit seqControl from above, since use same TTNA and trig
           %Event(n).seqControl = seq.seqControlIndTTNASD; % set ttna
        end
        
        n = n+1;
      end % separate B acquisition
      
      Event(n).info = 'Doppler acquisition.';
       if j==1 & evParms.state.SDOnlyInPlace
             %Event(n).tx = TXIndB; % experiment that did not work    
             Event(n).tx = TXIndSD;        
       else
             Event(n).tx = TXIndSD;        
       end
       
      Event(n).rcv = k+j+evParms.flag.separateBAcq;
      Event(n).recon = 0;      % no reconstruction.
      Event(n).process = 0;    % no processing
      Event(n).seqControl = seq.seqControlIndTTNASD; % set ttna
      if ismember(i, trig.frameNos) & ...
            ((j==1 & sendTriggers & ...
              trig.onFirstPulseOfTypeOnly & trig.onSD) | ...
             (sendTriggers & trig.onSD & ~trig.onFirstPulseOfTypeOnly))
        Event(n).seqControl = [Event(n).seqControl ...
                                 seq.seqControlIndTriggerOut];
      end
      n = n+1;
        
    end
    
    if lastFrameTTNA % try for multiangle B plane wave to see if
                     % affects frame rate
                     % now get the CFI frame
                     % replace last 2D acquisition Event's seqControl
      seqq = setdiff(Event(n-1).seqControl, seq.seqControlIndTTNASD);
      Event(n-1).seqControl = [seqq seq.seqControlIndTransferTTNA_postSD]; ...
          
    end
     
    if ~doCFI & syncFrameStart      
      Event(n-1).seqControl = 0;
    end
    
    % B/CFI events start here
    
    evStartCFI=n;

    %if ~singleTpcProfile
    %  Event(n-1).seqControl = [seq.seqControlIndSetTPCProfileCFI ...
         %                      seq.seqControlIndSetLongTTNAChangeProfile];
    %end

    % 20240117: shouldn't there be a long ttna after the last doppler?
    Event(n-1).seqControl = seq.seqControlIndTransferTTNA_postSD;

    Event(n).info = 'jump to end of CFI acqs';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no process
    %Event(n).seqControl = nsc; % test loopCnt
    Event(n).seqControl = seq.seqControlIndJumpOverCFIAcq(i);
    %SeqControl(nsc).command = 'jump'; 
    %seqControlIndPostCFIAcq = nsc; % need to set argument below
    %nsc=nsc+1;
    
    n=n+1;
    % this is point to jump to
    SeqControl(seq.seqControlIndLoopTstSkipCFI(i)).argument = n;

    if 0  % will crash VSX as standalone
    % always use this profile for doCFISeq
    Event(n).info = 'Set RcvProfile for B Doppler';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [seq.seqControlIndSetRcvProfile2D];
    n = n+1;
    end
    
    if evParms.flag.separateBAcq
      txCnt = 0;
      for j = 1:2 % two B-mode acqs
        txCnt=txCnt+1;
        Event(n).info = '2D acquisition.';
        if multiAngleB
          Event(n).tx = TXIndBVec(txCnt);  % use TX structure 1.
        else
          Event(n).tx = TXIndB;  % use TX structure 1.
        end
      
%        Event(n).rcv =  thisFrameRcvB(j); % use the two existing
                                          % rcvs for this frame:
                                          % bad idea
        Event(n).rcv = evParms.rcvOut.recvIndBCFIVec((i-1)*2+j);

        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
        Event(n).seqControl = seq.seqControlIndSetShortTTNACFI; % seq.seqControlIndTTNASD; % set ttna

        if j==1
          % set rcv profile before event
          Event(n).seqControl = [Event(n).seqControl];
          %seq.seqControlIndSetRcvProfile2D];
          %% this is kiss of death
        end
      
        if j==2
          Event(n).seqControl = [seq.seqControlIndSetTPCProfileCFI,...
                              seq.seqControlIndSetLongTTNAChangeProfile]; 
        end
        
        if (j==1 & sendTriggers & ...
            trig.onBPreCFI & trig.onFirstPulseofTypeOnly) | ...
              (sendTriggers & ...
               trig.onBPreCFI & ~trig.onFirstPulseofTypeOnly)
          
          Event(n).seqControl = [Event(n).seqControl ...
                              seq.seqControlIndTriggerOut];
        end
      
        n = n+1;
      end
      
    end
      
    Event(n).info = 'Set RcvProfile for CFI Doppler';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no processing
    Event(n).seqControl = [seq.seqControlIndSetRcvProfileCFI];
    n=n+1;

    % Acquire CFI Doppler ensemble.
    for j = 1:m % rays
      for k = 1:P(PIndCFI).dopPRIs
        Event(n).info = 'Acq CFI Doppler ensemble';
        Event(n).tx = TXOffset+j; % use next TX structure after 2D+SD.
        Event(n).rcv = ...
              rcvCFIStartInd(i)+(j-1)*P(PIndCFI).dopPRIs+k-1;
        Event(n).recon = 0;      % no reconstruction.
        Event(n).process = 0;    % no processing
                                 % TTNA for Doppler, sets intraensemble PRF
        Event(n).seqControl = seq.seqControlIndSetShortTTNACFI; 
        if ismember(i, trig.frameNos) & ((j==1 & k==1 & sendTriggers & ...
            trig.onFirstPulseOfTypeOnly & trig.onCFI) | ...
          (sendTriggers & trig.onCFI & ~trig.onFirstPulseOfTypeOnly))
          Event(n).seqControl = [Event(n).seqControl ...
                                 seq.seqControlIndTriggerOut];
        end
        n = n+1;
      end
    end

    if syncFrameStart
      % don't put a TTNA after last pulses if we are syncing at
      % start of frame using SW catch-up delay
      Event(n-1).seqControl = 0;
    end     

    if 1 | ~singleTpcProfile %| ~singleRcvProfile
                             % replace last Doppler acquisition 
                             % Event's seqControl
      seqq = setdiff(Event(n-1).seqControl, seq.seqControlIndSetShortTTNACFI);
      Event(n-1).seqControl = [seqq seq.seqControlIndSetTPCProfileSD,...
                               seq.seqControlIndTransferTTNA];
                          %seq.seqControlIndSetLongTTNAChangeProfile]; 
    end
    
     
    % set jump to after CFI sequential acqs if SD is running
    SeqControl(seq.seqControlIndJumpOverCFIAcq(i)).argument = n;
    
    % the two branches must share same tx2host
    Event(n).info = 'transfer acquired frame data to host';
    Event(n).tx = 0;         % no TX
    Event(n).rcv = 0;
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    %Event(n).seqControl = [nsc];
    Event(n).seqControl = seq.seqControlIndTransferToHost(i);
    %SeqControl(nsc).command = 'transferToHost';
    %nsc = nsc+1;
    n = n+1;
  
    % now do B recon for both branches
        
    Event(n).info = 'B recon';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = recon.reconIndB;  
    Event(n).process = 0;    % no process
    Event(n).seqControl = 0; % no seqCntrl

    if zeroBRecon
      Event(n).recon=0;
    end
    
    if evParms.flag.wtCapability 
      uscmsg(evParms.state.execState);
      n=n+1;   
      if evParms.flag.doWT
        % reconstruct PData for WT        
        Event(n).info = 'WTIQ recon';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = recon.reconIndWT;  
        Event(n).process = proc.procIndWT;    % no process
        Event(n).seqControl = 0; % no seqCntrl
      else
       
        % placeholder NOOP to keep sequence length constant between
        % SD/WT modes
        Event(n) = eventnoopfn(seq);
              
      end
      
    end
        
    if zeroBRecon
      Event(n).recon=0;
    end
         
    if evParms.flag.BLineWTCapability
      n=n+1;  
      % RF wall track processing
      if evParms.flag.BLineWT
        uscmsg(evParms.state.execState);
        Event(n).info = 'Process RF line';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;
        Event(n).process = proc.procIndBWT;
        Event(n).seqControl = 0; % no seqCntrl
      else
        Event(n) = eventnoopfn(seq);
      end    
    end
    
    n=n+1;
        
    % for CFI, jump over SD recon
    Event(n).info = 'loopTst jump over SD recon';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no process
%    Event(n).seqControl = nsc; % jump
    Event(n).seqControl = seq.seqControlIndJumpOverSDRecon(i);
    %SeqControl(nsc).command = 'loopTst';
    %SeqControl(nsc).condition = 'counter1';    
    %seqControlIndJumpOverSDRecon = nsc;
    %nsc=nsc+1;
    n=n+1;
    
    if ~ismember(i, evParms.largeSDParms.largeSDFrames)
      disp(['skipping frame ' num2str(i) ...
            ' for SDLarge recon at event ' num2str(n)] );
      lastSDRecon = evParms.gate.numSD;
      skipSDLarge=1;
       %lastSDRecon = evParms.gate.numSD+evParms.flag.largeSD;
    else
      skipSDLarge=0;
      lastSDRecon = evParms.gate.numSD+evParms.flag.largeSD;
    end
    
    % for B/SD mode, and B/SD/SDL
    Event(n).info = 'B image display'; 
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    
    if (evParms.flag.largeSD & ...
        evParms.largeSDParms.SDLargeDisplayOverlay)
      Event(n).process = proc.procIndB; % SDL will overlay later    
    else
      Event(n).process = proc.procIndBNow; % B only  
    end
    if zeroBRecon
      Event(n).process = 0;
    end
    
    
    Event(n).seqControl = 0; % no seqCntrl
    n=n+1;
    
    for q = 1:lastSDRecon
      Event(n).info = ['SD: recon and process ' num2str(q)];
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = recon.reconIndSDStart+q-1; % reconstruction
      if q <= evParms.gate.numSD | ...
         (evParms.flag.largeSD & ~ ...
          evParms.largeSDParms.doSDLargeInternalProc) 

      % if we are doing interal CDI, process in next event
        Event(n).process = proc.procIndSDExtStart+q-1;    % process
      else
        Event(n).process = 0;
      end      
      
      Event(n).seqControl = 0; % no seqCntrl
      n = n+1;
    end
    
    % kill sdlarge recon (MARK)
    %Event(n-1).recon = 0;
    
    if evParms.flag.largeSD & ...
       evParms.largeSDParms.doSDLargeInternalProc & ...
       ~skipSDLarge
      Event(n).info = 'Proc SDL-CDI';
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = 0;      % no reconstruction
      Event(n).process = proc.procIndSDLarge; 
      Event(n).seqControl = 0; % no seqCntrl
      n=n+1;
      
    end
    
    if evParms.flag.largeSD & evParms.largeSDParms.largeSDSave & ...
       ~skipSDLarge
      Event(n).info = 'Save SDL';
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = 0;      % no reconstruction
      Event(n).process = proc.procIndSDLargeSave; 
      Event(n).seqControl = 0; % no seqCntrl
      n=n+1;
    end
     
       
    if evParms.flag.largeSD & evParms.largeSDParms.procSDLargeOut & ...
       ~skipSDLarge
      Event(n).info = 'process CDI IQ output';
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = 0;      % no reconstruction
      if ~skipSDLarge
        Event(n).process = proc.procIndSDLargeOut;    % display
                                                      % processing
      else
        Event(n).process = proc.procIndSDLargeZeroOut;    
      end
      
      Event(n).seqControl = 0; % no seqCntrl
      n=n+1;
    end

    if evParms.flag.largeSD & ...
       evParms.largeSDParms.SDLargeDisplayOverlay
      Event(n).info = 'SDL overlay display';
      Event(n).tx = 0;         % no transmit
      Event(n).rcv = 0;        % no rcv
      Event(n).recon = 0;      % no reconstruction      
      Event(n).process = proc.procIndSDLargeIm;    % display processing
      Event(n).seqControl = 0; % no seqCntrl
      if isfield(evParms.flag, 'canSaveImage') & evParms.flag.canSaveImage
        n=n+1;
        Event(n).info = 'Image save';
        Event(n).tx = 0;         % no transmit
        Event(n).rcv = 0;        % no rcv
        Event(n).recon = 0;      % no reconstruction      
        Event(n).process = proc.procIndImSave;    % display processing
        Event(n).seqControl = 0; % no seqCntrl
      end
      
      n=n+1;
    end
     
    % for CFI, jump over SD recon
    Event(n).info = 'jump over CFI recon';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no process
    Event(n).seqControl =  seq.seqControlIndJumpOverCFIRecon(i); % jump
%    SeqControl(nsc).command = 'jump';
%    seqControlIndJumpOverCFIRecon = nsc;
%    nsc=nsc+1;
    
    n=n+1;
    
    Event(n).info = 'recon CFI';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = recon.reconIndCFI;  
    Event(n).process = 0;    
    Event(n).seqControl = 0;  
    SeqControl(seq.seqControlIndJumpOverSDRecon(i)).argument = n;
    n=n+1;
   
    Event(n).info = 'CFI Doppler processing';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = proc.procIndCFI; % IQ to estimate processing CFI
    Event(n).seqControl = 0; % no seqCntrl
    n = n+1;    

    Event(n).info = 'B image display';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = proc.procIndB;    % process 2D    
    Event(n).seqControl = 0; % no seqCntrl
    n=n+1;
    
    Event(n).info = 'CFI Doppler image display';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = proc.procIndCFIIm;    % display processing
    Event(n).seqControl = 0; % no seqCntrl
    n=n+1;

    % stay in CFI mode
    Event(n).info = 'set loopCnt back to 1';
    Event(n).tx = 0;         % no transmit
    Event(n).rcv = 0;        % no rcv
    Event(n).recon = 0;      % no reconstruction
    Event(n).process = 0;    % no process
    Event(n).seqControl = seq.seqControlIndSetLoopCnt1; % set loopCnt to 1    
    n=n+1;
     
    Event(n).info = 'exit to Matlab';
    Event(n).tx = 0;         % no TX
    Event(n).rcv = 0;
    Event(n).recon = 0;      % no reconstruction.
    Event(n).process = 0;    % no processing
    Event(n).seqControl = 0; % default
    if floor(i/3) == i/3     % Exit to Matlab every 3rd frame
        Event(n).seqControl = seq.seqControlIndReturnToMatlab;
    end
    SeqControl(seq.seqControlIndJumpOverCFIRecon(i)).argument = n;
    n = n+1;
    
end % frame 

Event(n).info = 'Jump back to first event';
Event(n).tx = 0;        % no TX
Event(n).rcv = 0;       % no Rcv
Event(n).recon = 0;     % no Recon
Event(n).process = 0;
Event(n).seqControl = seq.seqControlIndJumpToFirstEvent;
%SeqControl(nsc).command = 'jump';
SeqControl(seq.seqControlIndJumpToFirstEvent).argument = nStartAcq;

evParmsOut.ev.evIndStartSD = evIndStartSD;
evParmsOut.ev.nStartAcq = nStartAcq;
%evParmsOut.evIndFrameStart = evIndFrameStart;


end

