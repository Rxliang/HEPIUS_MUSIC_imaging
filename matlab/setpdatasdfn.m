function [PData, evParms] = setpdatasdfn(PData, x, z, evParms, changeParms)

if nargin < 5 | isempty(changeParms)
  changeParms.flag.changeRequested=0;
else
  changeParms.flag.changeRequested=1;
end

% x,z in wvl

gate = evParms.gate;
PDataIndB = gate.PDataIndB;
PDataIndSDStart = gate.PDataIndSDStart;
PDataIndSDEnd = gate.PDataIndSDEnd;
PDataIndCFI = gate.PDataIndCFI1;

% put first and last rays as close to edge as possible

transducerWidth_wvl = evParms.ev.numElementsUsed*evParms.Trans.spacing;

xAxSpacing_wvl = PData(PDataIndB).PDelta(1);
zAxSpacing_wvl = PData(PDataIndB).PDelta(3);

dopSDWidthRange_wvl = ...
    ceil(transducerWidth_wvl*evParms.FOV.SDRegionExtensionFactor/...
         xAxSpacing_wvl)* xAxSpacing_wvl;

imageWidth_wvl = evParms.P(evParms.ind.PIndB).imageWidth_mm / ...
                 evParms.P(evParms.ind.PIndB).lambda_mm;

imageWidth_mm = imageWidth_wvl * gate.lambda_mm;

if dopSDWidthRange_wvl > imageWidth_wvl 
  dopSDWidthRange_wvl = imageWidth_wvl;
  disp(['*** Warning: SD region allowed width larger than image width ' ...
        ' ***']);
end

zAxSpacing_wvl = PData(PDataIndB).PDelta(3);

% note, image width is not actually this
%ax.xAxBMode_wvl = -imageWidth_wvl/2:xAxSpacing_wvl:imageWidth_wvl/2;

xAxStart_wvl = PData(PDataIndB).Origin(1);
xAxEnd_wvl = xAxStart_wvl+xAxSpacing_wvl*(PData(PDataIndB).Size(2)-1);
ax.xAxBMode_wvl = xAxStart_wvl:xAxSpacing_wvl:xAxEnd_wvl;

% note, image height is not actually this
%ax.zAxBMode_wvl = evParms.FOV.dopR1:zAxSpacing_wvl:evParms.FOV.dopR2;

zAxStart_wvl = PData(PDataIndB).Origin(3);
zAxEnd_wvl = zAxStart_wvl+zAxSpacing_wvl*(PData(PDataIndB).Size(1)-1);
ax.zAxBMode_wvl = zAxStart_wvl:zAxSpacing_wvl:zAxEnd_wvl;

ax.xAxBMode_mm = ax.xAxBMode_wvl * gate.lambda_mm;
ax.zAxBMode_mm = ax.zAxBMode_wvl * gate.lambda_mm;

if ~evParms.gate.largeGateIQ 
  xSVWidth_wvl = gate.xSVWidth_mm/gate.lambda_mm;
  zSVWidth_wvl = gate.zSVWidth_mm/gate.lambda_mm;
  % actual used for spectrogram processing:  
  xSVWidthSpctr_wvl = xSVWidth_wvl; 
  zSVWidthSpctr_wvl = zSVWidth_wvl;   
else
  uscmsg(evParms.state.execState);
  xSVWidth_wvl = gate.xSVWidthIQ_mm/gate.lambda_mm;
  zSVWidth_wvl = gate.zSVWidthIQ_mm/gate.lambda_mm;
  % actual used for spectrogram processing:
  xSVWidthSpctr_wvl = gate.xSVWidth_mm/gate.lambda_mm;
  zSVWidthSpctr_wvl = gate.zSVWidth_mm/gate.lambda_mm;
end

% these are adjustments for gates placed beyond limits of image

indAdj = find(x-xSVWidth_wvl/2 < ax.xAxBMode_wvl(1));
x(indAdj) = ax.xAxBMode_wvl(1)+xSVWidth_wvl/2;

indAdj = find(x+xSVWidth_wvl/2> ax.xAxBMode_wvl(end));
x(indAdj) = ax.xAxBMode_wvl(end)-xSVWidth_wvl/2;

indAdj = find(z-zSVWidth_wvl/2 < ax.zAxBMode_wvl(1));
z(indAdj) = ax.zAxBMode_wvl(1)+zSVWidth_wvl/2;

indAdj = find(z+zSVWidth_wvl/2 > ax.zAxBMode_wvl(end));
z(indAdj) = ax.zAxBMode_wvl(end)-zSVWidth_wvl/2;

for PDataIndSD = PDataIndSDStart:PDataIndSDEnd
  zLen = max([1 round(zSVWidth_wvl/PData(PDataIndSD).PDelta(3))]);
  xLen = max([1 round(xSVWidth_wvl/PData(PDataIndSD).PDelta(1))]);

  ind = PDataIndSD-PDataIndSDStart+1;
  PData(PDataIndSD).Size(1,1) = zLen; % IQ samples for Doppler
  PData(PDataIndSD).Size(1,2) = xLen; % columns for Doppler
  PData(PDataIndSD).Size(1,3) = 1; % single image page

  xSDOriginOffset_wvl = x(ind) - xSVWidth_wvl/2;
  zSDOriginOffset_wvl = z(ind) - zSVWidth_wvl/2;

  disp(['PDataIndSD = ' num2str(PDataIndSD) ', Old origin: ' num2str(PData(PDataIndSD).Origin)]);
  PData(PDataIndSD).Origin = [xSDOriginOffset_wvl, 0, zSDOriginOffset_wvl];
  
  disp(['PDataIndSD = ' num2str(PDataIndSD) ', New origin: ' num2str(PData(PDataIndSD).Origin)]);
   
  gate.xSVWidth_wvl(ind) = xSVWidth_wvl;
  gate.zSVWidth_wvl(ind) = zSVWidth_wvl;

  gate.SD{ind}.xDelta_wvl = PData(PDataIndSD).PDelta(1);
  gate.SD{ind}.xEnd_wvl = xSDOriginOffset_wvl+...
                           gate.SD{ind}.xDelta_wvl*(xLen-1);   
  gate.SD{ind}.xAxSD_wvl = ...
      xSDOriginOffset_wvl:gate.SD{ind}.xDelta_wvl:gate.SD{ind}.xEnd_wvl;
  
  gate.SD{ind}.zDelta_wvl = PData(PDataIndSD).PDelta(3);  
  gate.SD{ind}.zEnd_wvl = zSDOriginOffset_wvl+...
                           gate.SD{ind}.zDelta_wvl*(zLen-1);   
  gate.SD{ind}.zAxSD_wvl = ...
      zSDOriginOffset_wvl:gate.SD{ind}.zDelta_wvl:gate.SD{ind}.zEnd_wvl;
  
  if evParms.gate.largeGateIQ 
    
    indCol = find((gate.SD{ind}.xAxSD_wvl >=  x(ind)-xSVWidthSpctr_wvl/2) ...
                  & (gate.SD{ind}.xAxSD_wvl <=  x(ind)+ ...
                     xSVWidthSpctr_wvl/2));
  
    indRow = find((gate.SD{ind}.zAxSD_wvl >=  z(ind)-zSVWidthSpctr_wvl/2) ...
                  & (gate.SD{ind}.zAxSD_wvl <=  z(ind)+ ...
                     zSVWidthSpctr_wvl/2));
    gate.SD{ind}.PDataSpectrRows = indRow;
    gate.SD{ind}.PDataSpectrCols = indCol;
    gate.lateralSubSampleIndex{ind} = indCol;
    gate.rangeSubSampleIndex{ind} = indRow;
  else
    gate.SD{ind}.PDataSpectrRows = 1:zLen;
    gate.SD{ind}.PDataSpectrCols = 1:xLen;
    gate.lateralSubSampleIndex{ind} = 1:xLen;
    gate.rangeSubSampleIndex{ind} = 1:zLen;

  end
 
end

% store these so we have them when we need to modify one at a time
gate.x = x;
gate.z = z;

PDataIndSDLarge = PDataIndSDEnd+1;
PDataIndSDLargeOut = PDataIndSDEnd+2;

gate.SDLarge.xSDDelta_wvl = 1;
gate.SDLarge.zSDDelta_wvl = 1;

if changeParms.flag.changeRequested & changeParms.flag.changeLargeSDWidth
  xWidth_mm = changeParms.changeLargeSDWidth_mm;
else
  xWidth_mm = gate.SDLarge.xSDWidth_mm;
end

dopSDWidthRange_mm = dopSDWidthRange_wvl*gate.lambda_mm;

if xWidth_mm > imageWidth_mm
  xWidth_mm = imageWidth_mm
  if changeParms.flag.changeRequested & changeParms.flag.changeLargeSDWidth
    disp(['*** warning: input width for CFI region larger than ' ...
          'maximum allowed. Setting to maximum. ***']);
  end
end

if changeParms.flag.changeRequested & ...
   changeParms.flag.changeLargeSDLateralOffset
  xOffset_mm = changeParms.changeLateralOffset_mm;
else
  xOffset_mm =  gate.SDLarge.xSDOffset_mm;
end

if changeParms.flag.changeRequested & ...
      changeParms.flag.changeLargeSDHeight
%  keyboard
  zHeight_mm = changeParms.largeSDHeight_mm;
else
  zHeight_mm = gate.SDLarge.zSDHeight_mm;
end

maxHeight_mm = ax.zAxBMode_mm(end)-ax.zAxBMode_mm(1);

if zHeight_mm > maxHeight_mm
  zHeight_mm = maxHeight_mm;
  if changeParms.flag.changeRequested & changeParms.flag.changeLargeSDHeight
    disp(['*** warning: input height for large SD region larger than ' ...
          'maximum allowed. Setting to maximum. ***']);
  end
end

if changeParms.flag.changeRequested & changeParms.flag.changeLargeSDZStart
  zStart_mm = changeParms.largeSDZStart_mm;
else
  zStart_mm = gate.SDLarge.zSDStart_mm;
end

if zStart_mm < ax.zAxBMode_mm(1)
  zStart_mm = ax.zAxBMode_mm(1);
  if changeParms.flag.changeRequested & changeParms.flag.changeLargeSDZStart
    disp(['*** warning: start depth of large SD region smaller than ' ...
          'minimum allowed. Setting to minimum. ***']);
  end
end

zEnd_mm= zStart_mm + zHeight_mm;

if zEnd_mm > ax.zAxBMode_mm(end)
  ax.zAxBMode_mm(end)
  zEnd_mm =  ax.zAxBMode_mm(end);
  zHeight_mm = zEnd_mm - zStart_mm 
  
  if changeParms.flag.changeRequested & changeParms.flag.changeLargeSDZStart
    disp(['*** warning: end depth of large SD region larger than ' ...
          'maximum allowed. Setting to maximum. ***']);    
  end  
end

zStartLim_mm = ax.zAxBMode_mm(end)-gate.SDLarge.zSDMinHeight_mm;
zEndLim_mm = ax.zAxBMode_mm(1)+gate.SDLarge.zSDMinHeight_mm;

if (zStart_mm > zStartLim_mm) | (zEnd_mm < zEndLim_mm)
  if zStart_mm > zStartLim_mm
    zStart_mm = zStartLim_mm;
  end

  if zEnd_mm < zEndLim_mm
    zEnd_mm = zEndLim_mm;
  end
  
  zHeight_mm =  zEnd_mm-zStart_mm;
else 
  zEnd_mm = zStart_mm +  zHeight_mm;  
end

% find B-image limits

xLeft_mm = xOffset_mm - xWidth_mm/2;
xRight_mm = xOffset_mm + xWidth_mm/2;

xLeftLim_mm = ax.xAxBMode_mm(end)-gate.SDLarge.xSDMinWidth_mm;
if xLeft_mm > xLeftLim_mm
  xLeft_mm = xLeftLim_mm;
end

xRightLim_mm = ax.xAxBMode_mm(1)+gate.SDLarge.xSDMinWidth_mm;
if xRight_mm < xRightLim_mm
  xRight_mm = xRightLim_mm;
end
 
if xLeft_mm < ax.xAxBMode_mm(1)
  disp(['*** warning: left side of largeSD region would be at ' ...
        num2str(xLeft_mm) ', left of ' ...
        'limit of ' num2str(ax.xAxBMode_mm(1))]);   
  xLeft_mm = ax.xAxBMode_mm(1);
end

xRight_mm = min(xLeft_mm + xWidth_mm,  ax.xAxBMode_mm(end));
  
if xRight_mm > ax.xAxBMode_mm(end)
  xRight_mm = ax.xAxBMode_mm(end);
  xLeft_mm = min(xRight_mm - xWidth_mm,  ax.xAxBMode_mm(1));
end

xWidth_mm = max(xRight_mm-xLeft_mm, gate.SDLarge.xSDMinWidth_mm);       
xOffset_mm = mean([xLeft_mm xRight_mm]);
  
gate.SDLarge.xSDWidth_mm = xWidth_mm;

gate.SDLarge.xSDOffset_mm = xOffset_mm;
gate.SDLarge.xSDStart_mm = gate.SDLarge.xSDOffset_mm - xWidth_mm/2;

xLeft_mm = gate.SDLarge.xSDStart_mm;
xRight_mm = gate.SDLarge.xSDStart_mm+xWidth_mm;

gate.SDLarge.zSDHeight_mm = zHeight_mm;
gate.SDLarge.zSDStart_mm = zStart_mm;
gate.SDLarge.xSDStart_wvl = xLeft_mm/gate.lambda_mm;
gate.SDLarge.xSDEnd_wvl = xRight_mm/gate.lambda_mm;

gate.SDLarge.zSDStart_wvl = zStart_mm/gate.lambda_mm;
gate.SDLarge.zSDEnd_wvl = zEnd_mm/gate.lambda_mm;

PData(PDataIndSDLarge).PDelta = [gate.SDLarge.xSDDelta_wvl, 0, ...
                                 gate.SDLarge.zSDDelta_wvl];
PData(PDataIndSDLarge).Origin = [gate.SDLarge.xSDStart_wvl, 0, ...
                                 gate.SDLarge.zSDStart_wvl];

ax.xAxSDL_wvl = gate.SDLarge.xSDStart_wvl:...
    gate.SDLarge.xSDDelta_wvl:gate.SDLarge.xSDEnd_wvl:...
    gate.SDLarge.xSDEnd_wvl;

ax.zAxSDL_wvl = gate.SDLarge.zSDStart_wvl:...
    gate.SDLarge.zSDDelta_wvl:gate.SDLarge.zSDEnd_wvl;

PData(PDataIndSDLarge).Size(1,1) = length(ax.zAxSDL_wvl); % IQ samples for Doppler
PData(PDataIndSDLarge).Size(1,2) = length(ax.xAxSDL_wvl); % columns for Doppler
PData(PDataIndSDLarge).Size(1,3) = 1; % single image page

PData(PDataIndSDLargeOut) = PData(PDataIndSDLarge);

ax.xAxSDL_mm = ax.xAxSDL_wvl*gate.lambda_mm;
ax.zAxSDL_mm = ax.zAxSDL_wvl*gate.lambda_mm;

[ax.XSDL_mm, ax.ZSDL_mm] = meshgrid(ax.xAxSDL_mm, ax.zAxSDL_mm);

gate.PDataIndSDLarge = PDataIndSDLarge;
gate.PDataIndSDLargeOut = PDataIndSDLargeOut;

gate.ax = ax;
evParms.gate = gate;
