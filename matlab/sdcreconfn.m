function [Recon, ReconInfo, recon, ri] = sdcreconfn(Recon, recon, ri, evParms, rcvParms)

doCFI = evParms.flag.doCFI;
numSD = evParms.gate.numSD;
multiAngleB = evParms.flag.multiAngleB;
numTXB = evParms.ev.numTXB;

if isempty(Recon)
  init=1;
else
  init=0;
end

if init
  recon = [];
  recon.reconIndB = 1;
  recon.reconIndCFI = 2;
  if doCFI % are we doing any form of CFI at all, at any time
    recon.reconIndSDStart = 3; %evParms.ind.reconIndSDStart;
    % B, spectral Doppler, CFI, largeSD, wall track
    numRecon = 1+numSD+1+evParms.flag.largeSD+evParms.flag.wtCapability; 
    recon.reconIndWT = numRecon;
  else
    recon.reconIndSDStart = 2;
    numRecon = 1+numSD; % B, spectral Doppler
  end

  if multiAngleB
    k = numTXB+1;
  else
    k = numTXB+2; % reconInfo(2) is for accumReplaceIQ
  end

  ri = [];
  if ~evParms.state.SDOnly % accommodate B mode ReconInfo 1st
    ri.RIStartIndSD(1) = k; 
  else
    ri.RIStartIndSD(1) = 1;
  end
  
  Recon = repmat(struct('senscutoff', evParms.reconParms.senscutoff, ...
                        'pdatanum', 1, ...
                        'rcvBufFrame', -1, ...  
                        'IntBufDest', [1,1], ...
                        'ImgBufDest', [1,-1], ...
                        'RINums', zeros(1,1)), 1, numRecon); 
  
end

if evParms.state.SDOnly | evParms.flag.doCFINow | ~ ...
      evParms.flag.separateBAcq | evParms.state.SDOnlyInPlace
  evParms.ev.priSkip=1;
end

numSDRecon =  numSD; 

for q = 1:numSDRecon
  ri.RIEndIndSD(q) = (ri.RIStartIndSD(q)+(evParms.ev.nPRIs/evParms.ev.priSkip)-1);
  if evParms.flag.largeSD | q < numSDRecon
    ri.RIStartIndSD(q+1) = ri.RIEndIndSD(q)+1;
  end
  
end


if evParms.flag.largeSD
  q=q+1;
  % to get fewer PRIs in recon for largeSD, need to divide total
  % PRIs (b+dop) by priSkip*skip for largeSD
  netSkip = (evParms.ev.priSkip*evParms.largeSDParms.priSkip);
  ri.RIEndIndSD(q) = (ri.RIStartIndSD(q)+...
                      floor(evParms.ev.nPRIs/netSkip)-1);

  %  ri.RIEndIndSD(q) = riLargeSDInd(evParms.largeSDParms.endPRI);

end


ri.RIEndInd = ri.RIEndIndSD(q); % advanced further if useCFI ccc
ri.SDIndVal = rcvParms.SDIndVal;    % used to allocate reconInfos
                                    % depending on SD mode in sdcreconinfofn

if doCFI
  m = evParms.P(rcvParms.PIndCFI).dopNumRays;
  ri.RIStartIndCFI = ri.RIEndInd+1;
  cfiPRIsRecond = evParms.P(rcvParms.PIndCFI).dopPRIs;
  ri.RIEndIndCFI = ri.RIStartIndCFI+m*cfiPRIsRecond-1;
end

% - Set Recon values for 2D frame.
if multiAngleB & init
  indRIUsedB = 1:numTXB;
  Recon(recon.reconIndB).RINums = [];
  Recon(recon.reconIndB).RINums(1,1:length(indRIUsedB)) = indRIUsedB;  % to increase
                                                                       % frame rate,
                                                                       % only recon
                                                                       % these
else
  Recon(recon.reconIndB).RINums(1,1:numTXB+1) = (1:numTXB+1);  % why do
                                                               % we need
                                                               % two,
                                                               % says
                                                               % need for
                                                               % synth
                                                               % aperture,
                                                               % but
                                                               % flashangles
                                                               % scripts
                                                               % do not
                                                               % have this
end

if doCFI 
  if init
    % Recon for CFI
    % - Set Recon values for Doppler ensemble.
    Recon(recon.reconIndCFI).pdatanum = evParms.ind.PDataIndCFI1; 
    Recon(recon.reconIndCFI).IntBufDest = [evParms.ind.interBufferIndCFI,1];
    Recon(recon.reconIndCFI).ImgBufDest = [0,0];
    %Recon(recon.reconIndCFI).senscutoff = 0.8;
  end
  % P.dopPRIs ReconInfos needed for Doppler ensemble.
  ri.numRICFI = m*cfiPRIsRecond;
  Recon(recon.reconIndCFI).RINums(1, 1:ri.numRICFI) = ...
      ri.RIStartIndCFI:ri.RIEndIndCFI;
  ri.RIEndInd=ri.RIEndIndCFI;
end

if 1 | ~evParms.state.SDOnly
  % first call to get ReconInfo unmodified by GUI 
  [ReconInfo, ri] = sdcreconinfofn(ri, evParms, rcvParms);
else
  evParms.ev.priSkip = 1;
  evParms.ev.numTXB=0;
  ri.RIEndIndSD =(ri.RIStartIndSD+(evParms.ev.nPRIs/evParms.ev.priSkip)-1);
  ri.RIEndInd = ri.RIEndIndSD(end); % advanced further if useCFI
                               % call to get new ReconInfo structure
  [ReconInfo,ri] = sdcreconinfofn(ri, evParms, rcvParms);
end

% - Set Recon values for spectral Doppler ensembles.

for i = 1:numSD 
  Recon(recon.reconIndSDStart+i-1).pdatanum = evParms.ind.PDataIndSDStart+i-1;
  Recon(recon.reconIndSDStart+i-1).IntBufDest = [evParms.ind.interBufferIndSDStart+i-1,1];
  Recon(recon.reconIndSDStart+i-1).ImgBufDest = [];

  % nPRIs/2 ReconInfos needed for Doppler.
  Recon(recon.reconIndSDStart+i-1).RINums = ...
      ri.RIStartIndSD(i):ri.RIEndIndSD(i);   
  Recon(recon.reconIndSDStart+i-1).senscutoff = 0.8; %reduce aperture to reduce angle spread
end

if evParms.flag.largeSD
  i=i+1;
  Recon(recon.reconIndSDStart+i-1).pdatanum = evParms.ind.PDataIndSDStart+i-1;
  Recon(recon.reconIndSDStart+i-1).IntBufDest = [evParms.ind.interBufferIndSDStart+i-1,1];
  Recon(recon.reconIndSDStart+i-1).ImgBufDest = [];
  Recon(recon.reconIndSDStart+i-1).RINums = ...
      ri.RIStartIndSD(i):ri.RIEndIndSD(i);   
  Recon(recon.reconIndSDStart+i-1).senscutoff = 0.8; 
end

if evParms.flag.wtCapability
  i=i+1;
  Recon(recon.reconIndWT).pdatanum = evParms.ind.PDataIndWT;
  Recon(recon.reconIndWT).IntBufDest = [evParms.ind.interBufferIndWT,1];
  Recon(recon.reconIndWT).ImgBufDest = [];
  Recon(recon.reconIndWT).RINums =  ri.RIStartIndWT:ri.RIEndIndWT;   
  Recon(recon.reconIndWT).senscutoff = 0.8; 
end


end

