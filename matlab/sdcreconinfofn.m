function [ReconInfo, ri] = sdcreconinfofn(ri, evParms, rcvParms)

doCFI = evParms.flag.doCFI;
doCFINow = 1; %evParms.flag.doCFINow;

TXIndSD = evParms.ev.TXIndSD;

% my understanding is separateBAcq means we actually acquire the B
% pulses and don't try to use SD to generate it.

if evParms.flag.separateBAcq
  TXIndB = evParms.ev.TXIndB;
else
  TXIndB = TXIndSD;
end

recvIndBVec = evParms.rcvOut.recvIndBVec;

% B-mode acq for the CFI sequence
recvIndBCFIVec = evParms.rcvOut.recvIndBCFIVec;
multiAngleB = evParms.flag.multiAngleB;
numTXB = evParms.ev.numTXB;

if evParms.flag.separateBAcq
  priSkip = evParms.ev.priSkip;
else
  priSkip = 1;
end

if ~evParms.state.SDOnly & evParms.flag.separateBAcq & ~evParms.state.SDOnlyInPlace
  priSkipSD = 2; % ReconInfos are generated for all modes, including
                 % duplex and CFI
else
  priSkipSD = 1;
end

if evParms.state.SDOnlyInPlace
  % preserve lower PRF for overlay 
  priSkipSDLargeReconInfo = 2;
else
  priSkipSDLargeReconInfo = priSkipSD*evParms.largeSDParms.priSkip;
end

priSkipRcv = evParms.ev.priSkipRcv;
nPRIs = evParms.ev.nPRIs;

if doCFINow
  TXIndCFI = evParms.ev.TXIndCFI;
  ri.RIStartIndCFI = ri.RIEndIndSD(end) + 1;  
  ri.RIEndIndCFI = ri.RIEndIndSD(end) + ri.numRICFI;
  ri.RIEndInd =  ri.RIEndIndCFI;
  % if WT-capable, will be increased below
end

% Define ReconInfo structures.
% - For 2D, we need 2 ReconInfo structures.
% - For Doppler, we only need nPRIs/2 ReconInfo structures for the
% default case, but we will define maxPRIs/2 structures to allow
% handling the maxPRI case.
% set replace IQ data as default:

ReconInfo = repmat(struct('mode', 'replaceIQ', ...    
    'txnum', 1, ...
    'rcvnum', 1, ...
    'pagenum',1, ...
    'regionnum', 1), [1 ri.RIEndInd]);

% - Set specific ReconInfo attributes.
%   - ReconInfos for 2D frame.

for k = 1:numTXB
  if k==1
    ReconInfo(k).mode = 'replaceIQ'; % replace IQ data
  else
    ReconInfo(k).mode = 'accumIQ'; % replace IQ data
  end
  
  ReconInfo(k).txnum = TXIndB+k-1;
  if ~evParms.state.SDOnly & ~evParms.state.SDOnlyInPlace
    % these Receives are not defined for SDOnly
    if evParms.flag.currentModeBCFI
        ReconInfo(k).rcvnum = recvIndBCFIVec(k);
    else
        ReconInfo(k).rcvnum = recvIndBVec(k);
    end
    
  end
  
end

if multiAngleB % not sure will work with CFI interleave mode
  ReconInfo(1).mode = 'replaceIQ';
  ReconInfo(numTXB).mode = 'accumIQ_replaceIntensity';

%  ReconInfo(numTXB+1).txnum = TXIndB+numTXB-1;  
%  ReconInfo(numTXB+1).rcvnum =  TXIndB+numTXB-1;
  k = numTXB; % offset of RINums after B mode
else
  k = numTXB+1; % offset of RINums after B mode
  ReconInfo(k).mode = 'accumIQ_replaceIntensity'; % accum and detect
  ReconInfo(k).txnum = TXIndB;  
  
  if evParms.flag.currentModeBCFI
    ReconInfo(k).rcvnum = recvIndBCFIVec(k);
  else
    ReconInfo(k).rcvnum = recvIndBVec(k);
  end

end

%if doCFINow % two adjacent Bs before rest
%  rcvSDOffset = ReconInfo(numTXB+1).rcvnum+1;
%else
 
if ~evParms.state.SDOnly & ~evParms.state.SDOnlyInPlace 
% 20240313 the above case appears flawed: when we have combo seq switchable
% to CFI mode, still want SD reconinfos to contain receives at sequence start
    % when we have interleved b and sd seqs, we need to
    % start one before 2nd B-mode recv
    %    rcvSDOffset = ReconInfo(numTXB+1).rcvnum-1;
    rcvSDOffset = 2;
else
  rcvSDOffset = 2;
end

if evParms.state.SDOnlyInPlace 
  rcvSDOffset=2;
end

% this appears to be equal to the total number of SD recons
for q = 1:evParms.gate.numSD
  k = ri.RIStartIndSD(q)-1; 
  % why do we have this, this is for B:
  %ReconInfo(k).mode = 'accumIQ_replaceIntensity';
  %  - ReconInfos for SD Doppler ensemble.
  % changed 20240302: ro priSkipRcv, which should be 2 for interleaved b
  for j = 1:nPRIs/priSkip
      ReconInfo(k+j).mode = 'replaceIQ';
    if (j==1 | j==2) & evParms.state.SDOnlyInPlace
      ReconInfo(k+j).txnum = TXIndB;
      ReconInfo(k+j).txnum = TXIndSD; % experiment that did not work      
    else
      ReconInfo(k+j).txnum = TXIndSD;
    end
    ReconInfo(k+j).rcvnum = priSkipSD*(j-1)+rcvSDOffset- ...
        (evParms.state.SDOnly | evParms.state.SDOnlyInPlace); 
    %ReconInfo(k+j).rcvnum = priSkipSD*(j-1)+rcvSDOffset;    
    ReconInfo(k+j).pagenum = j;
    % nPRIs/2 pages in SD interbuffer
    ReconInfo(k+j).regionnum = 1;  % only 1 col..
  end
end

if evParms.flag.largeSD
   q = q+1;
   k = ri.RIStartIndSD(q)-1; 
   J = ri.RIEndIndSD(q)-ri.RIStartIndSD(q)+1;
   for j = 1:J
    ReconInfo(k+j).mode = 'replaceIQ';
    ReconInfo(k+j).txnum = TXIndSD;
    if evParms.flag.separateBAcq
       ReconInfo(k+j).rcvnum = priSkipSDLargeReconInfo*(j-1)+rcvSDOffset- ...
           evParms.state.SDOnly; 
    else
       ReconInfo(k+j).rcvnum = evParms.largeSDParms.priSkip*(j-1)+ ...
           rcvSDOffset;    
    end
    
    % destination page in interbuffer
    ReconInfo(k+j).pagenum = j;
    % nPRIs/2 pages in SD interbuffer
    ReconInfo(k+j).regionnum = 1;  % only 1 col..
  end
end

if evParms.flag.wtCapability
  ri.RIStartIndWT = ri.RIEndIndCFI+1;
  numWTTx = 1;
  ri.RIEndIndWT = ri.RIStartIndWT+numWTTx-1;  
  J=numWTTx;
  k=ri.RIStartIndWT;
   for j = 1:J
    ReconInfo(k+j-1).mode = 'replaceIQ';
    ReconInfo(k+j-1).txnum = TXIndB;
    if evParms.flag.separateBAcq
     if ~evParms.state.SDOnlyInPlace
       ReconInfo(k+j-1).rcvnum =  recvIndBVec(1);       
     else
       ReconInfo(k+j-1).rcvnum =  1; % placeholder, since no
                                     % B-receives in SDOnlyInPlace
                                     % mode       
     end
     
    end    
    % destination page in interbuffer
    ReconInfo(k+j-1).pagenum = j;
    % nPRIs/2 pages in SD interbuffer
    ReconInfo(k+j-1).regionnum = 1;  % only 1 col..
  end
end

% now add in CFI ReconInfos
% CFI: Define ReconInfo structures.
%  - ReconInfos for CFI Doppler ensemble.

if doCFINow
  ReconInfo(ri.RIStartIndCFI).Pre = 'clearInterBuf';
  TXOffset=TXIndCFI-1;; % one for all B-mode, next for all spectral
                        % Doppler, next for CFI

  if ~evParms.flag.doCFISeq 
    rcvCFIInd = find(evParms.rcvOut.rcvTypeInd==rcvParms.CFIIndVal);
    rcvCFIStartInd =  rcvCFIInd(1); % just used for TOF calc
  end
  
  m = evParms.P(rcvParms.PIndCFI).dopNumRays;
  cfiPRIsRecond = evParms.P(rcvParms.PIndCFI).dopPRIs;
    for n = 1:m 
      offset = ri.RIStartIndCFI-1+n;
      for j = 1:cfiPRIsRecond
        ReconInfo(offset+(j-1)*m).threadSync = 1;
        ReconInfo(offset+(j-1)*m).mode = 'replaceIQ';
        ReconInfo(offset+(j-1)*m).txnum = TXOffset + n;
%        if ~evParms.flag.doCFISeq
%          ReconInfo(offset+(j-1)*m).rcvnum = rcvCFIInd(j);
%        else
          % index here is per frame
          ReconInfo(offset+(j-1)*m).rcvnum = rcvCFIStartInd(1)+(n-1)*cfiPRIsRecond+j-1;
%        end
      
        ReconInfo(offset+(j-1)*m).regionnum = n; % nth Doppler ray
        ReconInfo(offset+(j-1)*m).pagenum = j;
        ReconInfo(offset+(j-1)*m).scaleFactor = 0.2;
       % if n==2
       % keyboard
       % end
        
    end
  end
end


end

