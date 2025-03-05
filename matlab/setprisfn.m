function [P, nPRIsReduced, nDopPRIsUsed] = setprisfn(evParms)

doCFI = evParms.flag.doCFI;
doCFISeq = evParms.flag.doCFISeq;
doCFISep = evParms.flag.doCFISep;
PIndCFI = evParms.ind.PIndCFI;
nPRIs = evParms.ev.nPRIs;
P = evParms.P;

% need to keep number of total PRIs the same when adding CFI, so
% need to reduce B/SD pairs
if doCFISeq
  % reduce B/SD pairs by number of CFI pulses
  %nPRIsReduced = nPRIs -  m*P(PIndCFI).dopPRIs; % this is if we
  %can reduce the B/SD to include sequential CFI, but this has
  %issue because of CFI recon time
  nPRIsReduced = nPRIs;
  if ~iseven(nPRIsReduced)
    error('nPRIsReduced is not even.');
  end
else
  nPRIsReduced = nPRIs;
end

% spectralDoppler function needs this to get number of pulses
if doCFI & ~doCFISeq & ~doCFISep
  nDopPRIsUsed = nPRIsReduced/3;
  P(PIndCFI).dopPRIs = nDopPRIsUsed;
else
  if evParms.flag.separateBAcq
    nDopPRIsUsed = nPRIsReduced/2;
  else
    nDopPRIsUsed = nPRIsReduced;
  end
end
