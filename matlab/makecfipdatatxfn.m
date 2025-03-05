function [PData, TX, evParms] = makecfipdatatxfn(PData, TX, evParms)

% make PData and TX structures for CFI

PDataIndB = evParms.ind.PDataIndB;
PDataIndCFI1 = evParms.ind.PDataIndCFI1;
PDataIndCFI2 = evParms.ind.PDataIndCFI2;
P = evParms.P;
PIndB = evParms.ind.PIndB;
PIndCFI = evParms.ind.PIndCFI;

% put first and last rays as close to edge as possible
firstRayLocX = -(ceil( (P(PIndCFI).dopNumTx-1)/2))*evParms.Trans.spacing;
lastRayLocX = ceil((P(PIndCFI).dopNumTx-1)/2)*evParms.Trans.spacing;
TxDopOrgX = linspace(firstRayLocX, lastRayLocX, P(PIndCFI).dopNumRays);
minEl = min(evParms.ind.txUsedInd);
maxEl = max(evParms.ind.txUsedInd);

[dum, elInd] = ...
    findclosestinvec(evParms.Trans.ElementPos(evParms.ind.txUsedInd,1), ...
                     TxDopOrgX);
P(PIndCFI).dopTxOrgChnl = evParms.ind.txUsedInd(elInd);
imageWidth = P(evParms.ind.PIndB).imageWidth_mm / P(evParms.ind.PIndB).lambda_mm;
transducerWidth = length(evParms.ind.txUsedInd)*evParms.Trans.spacing;
dopPDataWidth = transducerWidth*P(PIndCFI).regionExtensionFactor;

if P(PIndCFI).sameShapeAsB
  PData(PDataIndCFI1) = PData(PDataIndB);
else
  % set origin the same as B mode, then r1 and r2 with index from there
  PData(PDataIndCFI1).Origin = PData(PDataIndB).Origin;
  PData(PDataIndCFI1).Origin(1) = -dopPDataWidth/2;
  PData(PDataIndCFI1).Origin(3) = P(PIndCFI).dopStartDepth;
end

PData(PDataIndCFI1).PDelta(1) = evParms.Trans.spacing/2;
PData(PDataIndCFI1).PDelta(2) = PData(PDataIndB).PDelta(2);
PData(PDataIndCFI1).PDelta(3) = PData(PDataIndB).PDelta(3);

PData(PDataIndCFI1).Size(1) = ...
    ceil((P(PIndCFI).dopEndDepth-P(PIndCFI).dopStartDepth)/PData(PDataIndCFI1).PDelta(3));
PData(PDataIndCFI1).Size(2) = ceil(dopPDataWidth ...
                                    / PData(PDataIndCFI1).PDelta(1)); 
PData(PDataIndCFI1).Size(3) = PData(PDataIndB).Size(3);

m = P(PIndCFI).dopNumRays;

%parallelogramWid = P(PIndCFI).regionOverlapFactor*...
%    P(PIndCFI).dopDispEle*evParms.Trans.spacing/m;

parallelogramWid = (dopPDataWidth/m)*P(PIndCFI).regionOverlapFactor;

%parallelogramWid = P(PIndCFI).dopNumTx*evParms.Trans.spacing/m;

PData(PDataIndCFI1).Region = repmat(struct('Shape',struct(...
    'Name','Parallelogram',...
    'Position',[0,0,P(PIndCFI).dopStartDepth],... %wil be changed later
    'width', parallelogramWid,...
    'height', P(PIndCFI).dopEndDepth-P(PIndCFI).dopStartDepth,...
    'angle',P(PIndCFI).dopAngle)),1, m);

for n = 1:m    
    PData(PDataIndCFI1).Region(n).Shape.Position(1) = TxDopOrgX(n);
end

% expand edge regions to edge of PData(PDataIndCFI1). What is new midpoint for
% parallelogram?

expandCFIToEdges=P(PIndCFI).expandCFIToEdges.flag;
if expandCFIToEdges
  REdge1 = TxDopOrgX(1)+parallelogramWid/2;
  wt = P(PIndCFI).expandCFIToEdges.weight;
  wDiff = REdge1-PData(PDataIndCFI1).Origin(1);
  newLPGramCenter = PData(PDataIndCFI1).Origin(1)+ wt * wDiff; 
  LWid = (REdge1-newLPGramCenter)*2;
 
  PData(PDataIndCFI1).Region(1).Shape.Position(1) = newLPGramCenter;
  PData(PDataIndCFI1).Region(1).Shape.width = LWid;
  
  LEdge2 = TxDopOrgX(end)-parallelogramWid/2;
  xEnd = PData(PDataIndCFI1).Origin(1)+PData(PDataIndCFI1).Size(2)*PData(PDataIndCFI1).PDelta(1);
  wDiff = xEnd-LEdge2;
  newRPGramCenter = LEdge2 +(1-wt) * wDiff;
  RWid = (newRPGramCenter-LEdge2)*2;
  PData(PDataIndCFI1).Region(n).Shape.Position(1) = newRPGramCenter;
  PData(PDataIndCFI1).Region(n).Shape.width = RWid;
  
end

PData(PDataIndCFI1).Region = computeRegions(PData(PDataIndCFI1));

% make region enclosing only CFI PData regions, assumes dopAngle==0

regEdgeL = [];
regEdgeR = [];
for n = 1:m
  regEdgeL(n) = PData(PDataIndCFI1).Region(n).Shape.Position(1)- ...
      PData(PDataIndCFI1).Region(n).Shape.width/2;
  regEdgeR(n) = PData(PDataIndCFI1).Region(n).Shape.Position(1)+ ...
      PData(PDataIndCFI1).Region(n).Shape.width/2;  
end
regEdge = [regEdgeL regEdgeR];
regWidCFI = max(regEdge)-min(regEdge);
evParms.P(PDataIndCFI1).regWidCFI=regWidCFI;
evParms.P(PDataIndCFI1).maxLateralExtentCFI=max(abs(regEdge)); % for
                                                               % receive
                                                               % length calc

% - PData(PDataIndCFI2) is used to outline the doppler area
PData(PDataIndCFI2) = PData(PDataIndCFI1); PData(PDataIndCFI2).Region = [];
% this seems wrong, shouldn't this have origin(1) same as CFI1?
PData(PDataIndCFI2).Region.Shape = struct(...
    'Name','Parallelogram',...
    'Position',[PData(PDataIndCFI2).Origin(1)+...
                PData(PDataIndCFI2).Size(2)*PData(PDataIndCFI2).PDelta(1)/2,0,...
                P(PIndCFI).dopStartDepth],...
    'width', regWidCFI,...
    'height', P(PIndCFI).dopEndDepth-P(PIndCFI).dopStartDepth,...
    'angle',P(PIndCFI).dopAngle);
PData(PDataIndCFI2).Region = computeRegions(PData(PDataIndCFI2));

winNum = P(PIndCFI).dopNumTx*2;
W = hannfn(winNum)'; % 128 elements for hann window
ApodRange = round(-P(PIndCFI).dopNumTx/2:P(PIndCFI).dopNumTx/2-1);

TXCFIOffset = P(PIndB).numRays;

numDelays = evParms.Trans.numelements;
apodNotSet =  zeros(1, numDelays);
apodSet =  ones(1, numDelays);

if evParms.flag.BLineWT
  uscmsg(evParms.state.execState);
  numElTxBWT = evParms.BLineWT.numElTxBWT;
  TxBWTOrgX = [];
  ApodRangeBWT = [];
  for i = 1:length(evParms.BLineWT.originTX_wvl)
    TxBWTOrgX(i) = evParms.BLineWT.originTX_wvl;
    [dum, TxBWTOrgChnl] = ...
      findclosestinvec(evParms.Trans.ElementPos(:,1), ...
                       TxBWTOrgX(i));
  end
  ApodRangeBWT = round(-numElTxBWT/2:numElTxBWT/2-1);
  numTx =  TXCFIOffset+m+evParms.BLineWT.numOrigin;
  
else
  numTx =  TXCFIOffset+m;
end

% Specify TX structure array.  
if evParms.flag.useAperture
  TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', apodNotSet, ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,numDelays),...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'aperture', 1, ...                   
                   'peakBLMax', 15.0), 1, numTx);
else
  TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', apodNotSet, ...
                   'focus', 0, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,numDelays),...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax', 15.0), 1, numTx);
end

TW = evParms.TW;

for n = 1:m
    TX(TXCFIOffset+n).Origin(1) = TxDopOrgX(n);  
    % Apodize existing channels only
    ApodPreInd = P(PIndCFI).dopTxOrgChnl(n)+ApodRange;
    ind = find(ApodPreInd >= minEl & ApodPreInd <= ...
               maxEl);
    winPre = W(winNum/2+ApodRange);
    if ~evParms.flag.zeroTxApod
      TX(TXCFIOffset+n).Apod(ApodPreInd(ind)) = winPre(ind);
    else
      TX(TXCFIOffset+n).Apod(ApodPreInd(ind)) = 0;        
    end
    
    TX(TXCFIOffset+n).waveform = evParms.ind.TWIndCFI;
    TX(TXCFIOffset+n).Steer = [P(PIndCFI).dopAngle,0.0];
    TX(TXCFIOffset+n).focus = P(PIndCFI).txFocus;
    TX(TXCFIOffset+n).Delay = computeTXDelays(TX(TXCFIOffset+n));
    % computeTXPD reads global TW, or crashes. Yuck.
    TX(TXCFIOffset+n).TXPD = computeTXPD(TX(TXCFIOffset+n), ...
                                         PData(PDataIndCFI1));

end

if evParms.flag.BLineWT
  uscmsg(evParms.state.execState);
  TXBWTOffset = TXCFIOffset+m;
  for n = 1:evParms.BLineWT.numOrigin
    TX(TXBWTOffset+n).Origin(1) = TxBWTOrgX(n);  
    % Apodize existing channels only
    ApodPreInd = TxBWTOrgChnl(n)+ApodRangeBWT;
    %ind = find(ApodPreInd >= min( TxBWTOrgChnl) & ApodPreInd <= ...
    %           max(TxBWTOrgChnl));
    ind = find(ApodPreInd >= 1 & ApodPreInd <= numDelays);    
    winPre = W(winNum/2+ApodRangeBWT);
    if ~evParms.flag.zeroTxApod
      TX(TXBWTOffset+n).Apod(ApodPreInd(ind)) = winPre(ind);
    else
      TX(TXBWTOffset+n).Apod(ApodPreInd(ind)) = 0;
    end
            
    TX(TXBWTOffset+n).waveform = evParms.ind.TWIndBWT;
    TX(TXBWTOffset+n).Steer = [0,0.0];
    TX(TXBWTOffset+n).focus = evParms.BLineWT.focus_wvl;
%        keyboard
    TX(TXBWTOffset+n).Delay = computeTXDelays(TX(TXBWTOffset+n));
    % computeTXPD reads global TW, or crashes. Yuck.
    %TX(TXBWTOffset+n).TXPD = computeTXPD(TX(TXBWTOffset+n), ...
    %                                     PData(PDataIndCFI1));
  end
  evParms.ev.TXIndBWTStart = TXBWTOffset+1;
end



%keyboard

end

