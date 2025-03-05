function [Receive, rcvOut] = sdcallreceivefn(rcvParms, evParms, P, Receive)

if nargin < 4
  Receive = [];
end

% All receives needed are allocated once and for all
% When switching to SDOnly, the relevant receives are modified

% think it better to reallocate them
allocateMax = 0; 

if evParms.state.SDOnly
  evParms.ev.priSkip = 1;
end

nfrms = evParms.ev.nfrms;
priSkip = evParms.ev.priSkip;
maxPRIs = evParms.ev.maxPRIs;
nPRIs = evParms.ev.nPRIs;

if ~iseven(maxPRIs)
  error('maxPRIs needs to be even so we have B and SD pairs');
end

m = P(rcvParms.PIndCFI).dopNumRays;

if evParms.state.SDOnly
  numRcvBSDPerFrame = nPRIs; % allocate only as needed
  numRcvPerFrame = numRcvBSDPerFrame;
else
  if allocateMax
    numRcvBSDPerFrame = maxPRIs; % allocate maximum. okay to use < max
                               % Receive in a sequence
  else
    numRcvBSDPerFrame = nPRIs;
  end
  numRcvBCFIPerFrame = 2;
  numRcvCFIPerFrame = m*P(rcvParms.PIndCFI).dopPRIs+numRcvBCFIPerFrame;
  numRcvPerFrame = numRcvBSDPerFrame+numRcvCFIPerFrame;
end

Receive = [];

numRcv = nfrms*numRcvPerFrame;

if isempty(Receive) | evParms.state.SDOnly
  if evParms.flag.useAperture
    Receive = repmat(struct('Apod', rcvParms.apodNotSet, ...
                        'startDepth', P(rcvParms.PIndB).startDepth, ...
                        'endDepth', 0, ...
                        'TGC', rcvParms.TGCIndB, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'aperture', 1, ... % MUX
                        'sampleMode',  rcvParms.sampleMode.B,...
                        'mode', 0, ...
                        'callMediaFunc', 0),1, numRcv);
  else
    Receive = repmat(struct('Apod', rcvParms.apodNotSet, ...
                        'startDepth', P(rcvParms.PIndB).startDepth, ...
                        'endDepth', 0, ...
                        'TGC', rcvParms.TGCIndB, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode',  rcvParms.sampleMode.B,...
                        'mode', 0, ...
                        'callMediaFunc', 0),1, numRcv);
  end

end
   
% start of each CFI ensemble
rcvCFIStartInd = [];
rcvFrameStartInd = [];
rcvInd = 0;
rcvBCnt = 0;
rcvBCFICnt = 0; % B's before CFI
rcvOut.recvIndBVec = [];
maxLateralDop = evParms.P(evParms.ind.PDataIndCFI1).maxLateralExtentCFI;

if ~evParms.flag.separateBAcq | evParms.state.SDOnlyInPlace
  priSkip=1;
end

%priSkip=2; % always account for B/SD

if ~evParms.state.SDOnly
  BSDPrisPerFrame = maxPRIs;
else
  BSDPrisPerFrame = nPRIs;
end

if ~allocateMax
  BSDPrisPerFrame = nPRIs;
end
  
rcvInd=1;
for i = 1:nfrms
  rcvFrameStartInd(i) = rcvInd;
  %rcvFrameStartInd(1:nfrms) = rcvInd+1;
  k = rcvInd-1;
  acqCnt=0;
  for j = 1:priSkip:BSDPrisPerFrame
    % keep Receives same for SDOnly or not
    if ~evParms.state.SDOnly & evParms.flag.separateBAcq & ~evParms.state.SDOnlyInPlace
      % -- 2D acquisition.
      Receive(k+j).Apod(rcvParms.indRecv) = 1.0;
      Receive(k+j).framenum = i;
      Receive(k+j).acqNum = j;
      Receive(k+j).sampleMode =  rcvParms.sampleMode.B;
      % need to set this, because VSX can adjust these, and each
      % Recon needs all receives to have same endDepth
      Receive(k+j).startDepth = P(rcvParms.PIndB).startDepth;      
      Receive(k+j).endDepth = evParms.FOV_B.maxDist_wvl;
      Receive(k+j).demodFrequency = rcvParms.TW(evParms.ind.TWIndB).Parameters(1);
      rcvTypeInd(rcvInd)=rcvParms.BIndVal;
      rcvBCnt = rcvBCnt+1;
      rcvOut.recvIndBVec(rcvBCnt) = rcvInd;
      rcvInd = rcvInd+1;
      
      if evParms.flag.BLineWT && ismember(j, evParms.BLineWT.TXIndBWT)
        Receive(k+j).demodFrequency = rcvParms.TW(evParms.ind.TWIndBWT).Parameters(1);
        %        [j k+j]
      end
      
      if 0 & evParms.flag.BLineWT & j~=1 & j ~=2 & j~=3
          uscmsg(evParms.state.execState);
          rayIndBWT = evParms.BLineWT.numTxBWT - (j - 4);
          if rayIndBWT >=0
             Receive(k+j).sampleMode = 'custom';
             Receive(k+j).ADCRate = 35.7143;             
          end
      end
       
      if j==1, Receive(k+j).callMediaFunc = 1; end
      % -- Doppler.
      % start Dop aperture at element 33
      Receive(k+j+1).Apod(rcvParms.indRecvDoppler) = 1.0; 
      Receive(k+j+1).TGC = rcvParms.TGCIndSD;
      Receive(k+j+1).sampleMode = rcvParms.sampleMode.SD;
      Receive(k+j+1).startDepth = P(rcvParms.PIndSD).startDepth;      
      Receive(k+j+1).endDepth = evParms.FOV_SD.maxDist_wvl;
      Receive(k+j+1).framenum = i;
      if evParms.state.SDOnly
        acqCnt = acqCnt+1;
        Receive(k+j+1).acqNum = acqCnt;
        Receive(k+j).acqNum = acqCnt; % make sure have sequential order
      else
        Receive(k+j+1).acqNum = j+1;
      end
      rcvTypeInd(rcvInd)=rcvParms.SDIndVal;      
      rcvInd = rcvInd+1;
      acqOffset = j+1; % because last Receive was k+j+1    
      
    else % SDOnly (also SDOnlyInPlace). No B mode at all for main sequence
       if j==1, Receive(k+j).callMediaFunc = 1; end
       % -- Doppler.
       Receive(k+j).Apod(rcvParms.indRecvDoppler) = 1.0; 
       Receive(k+j).TGC = rcvParms.TGCIndSD;
       Receive(k+j).sampleMode = 'BS100BW';
       Receive(k+j).startDepth = P(rcvParms.PIndSD).startDepth;      
       Receive(k+j).endDepth = evParms.FOV_SD.maxDist_wvl;
       Receive(k+j).framenum = i;
       Receive(k+j).acqNum = j;
       % need to set this, because VSX can adjust these, and each
       % Recon needs all receives to have same endDepth
       rcvTypeInd(rcvInd)=rcvParms.SDIndVal;
       rcvInd = rcvInd+1;
       if ~evParms.flag.separateBAcq % B needs to be reconstructed
         rcvBCnt = rcvBCnt + 1;
         rcvOut.recvIndBVec(rcvBCnt) = rcvInd;
       end
       acqOffset = j;    
    end % SDOnly

  end % all B/SD PRIs

%  acqOffset = j+evParms.flag.separateBAcq; % no 2nd receive unless
                                           % there is B

  
  % not needed to generate these again, even if SDOnly
  if ~evParms.state.SDOnly
    
     % need two B mode receives for B-mode that precedes CFI
     for q = 1:2
       Receive(rcvInd).Apod(rcvParms.indRecv) = 1.0;
       Receive(rcvInd).framenum = i;
       Receive(rcvInd).acqNum = acqOffset+q;
       Receive(rcvInd).sampleMode =  rcvParms.sampleMode.B;
       Receive(rcvInd).startDepth = P(rcvParms.PIndB).startDepth;      
       Receive(rcvInd).endDepth = evParms.FOV_B.maxDist_wvl;
       Receive(rcvInd).demodFrequency = rcvParms.TW(evParms.ind.TWIndB).Parameters(1);
       rcvTypeInd(rcvInd)=rcvParms.BIndVal;
       rcvBCFICnt = rcvBCFICnt +1;
       rcvOut.recvIndBCFIVec(rcvBCFICnt) = rcvInd;
       rcvInd = rcvInd+1;
     end
      
    rcvCFIStartInd(i) = rcvInd;
    acqOffset = acqOffset+q;
    
    for r = 1:m
      for k = 1:P(rcvParms.PIndCFI).dopPRIs
        Receive(rcvInd).Apod(rcvParms.indRecv) = 1.0;
        Receive(rcvInd).framenum = i;
        if evParms.state.SDOnly
          Receive(rcvInd).acqNum = acqCnt+(r-1)* ...
              P(rcvParms.PIndCFI).dopPRIs+k;
        else
          Receive(rcvInd).acqNum = acqOffset+(r-1)* ...
              P(rcvParms.PIndCFI).dopPRIs+k;
        end
        Receive(rcvInd).TGC = rcvParms.TGCIndCFI;    
        Receive(rcvInd).startDepth = P(rcvParms.PIndCFI).startDepth;
        Receive(rcvInd).endDepth = evParms.FOV_CFI.maxDist_wvl;
        Receive(rcvInd).InputFilter = rcvParms.BPF;
        Receive(rcvInd).demodFrequency = rcvParms.TW(rcvParms.TWIndCFI).Parameters(1);
        Receive(rcvInd).sampleMode = 'BS100BW';
        rcvTypeInd(rcvInd)=rcvParms.CFIIndVal;
        rcvInd=rcvInd+1;
      end
    end
  end
  
end % frames

rcvOut.rcvTypeInd = rcvTypeInd;
rcvOut.rcvFrameStartInd = rcvFrameStartInd;
rcvOut.rcvCFIStartInd = rcvCFIStartInd;
rcvOut.numRcv = numRcv;

% size check
if 0
  size([Receive.framenum])
  size([Receive.acqNum])
  size([Receive.endDepth])
  size([Receive.bufnum])
  size([Receive.mode])
end



