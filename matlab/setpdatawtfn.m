function [PData, evParms] = setpdatawtfn(PData, x, z, evParms, changeParms)

if nargin < 5 | isempty(changeParms)
  changeParms.flag.changeRequested=0;
else
  changeParms.flag.changeRequested=1;
end

disp('in setpdatawtfn')
% x,z in wvl

gate = evParms.gate;
PDataIndB = gate.PDataIndB;
PDataIndWT = evParms.ind.PDataIndWT;

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

% if ~evParms.gate.largeGateIQ 
if 1
  %uscmsg(evParms.state.execState);
%  xSVWidth_wvl = gate.xSVWidthIQ_mm/gate.lambda_mm;
%  zSVWidth_wvl = gate.zSVWidthIQ_mm/gate.lambda_mm;
  
   xSVWidthPre_wvl = evParms.wt.xWidth_mm/gate.lambda_mm;
   zSVWidthPre_wvl = evParms.wt.zHeight_mm/gate.lambda_mm;
  
end

% this code sets PData based on x,z marker, rather than menu.
% below is the ui-adjustable version

PData(PDataIndWT) = [];


if 1
  ind = 1; % for now fix on (instance)
  xWTOriginOffset_wvl = x(ind) - xSVWidthPre_wvl/2;
  % restrict to B-mode limits
  xWTOriginOffset_wvl = max(xWTOriginOffset_wvl, xAxStart_wvl);
  
  zWTOriginOffset_wvl = z(ind) - zSVWidthPre_wvl/2;
  zWTOriginOffset_wvl = max(zWTOriginOffset_wvl, zAxStart_wvl);
  
  xWidthMax_wvl = xAxEnd_wvl-xWTOriginOffset_wvl;
  xSVWidth_wvl = min(xSVWidthPre_wvl,  xWidthMax_wvl);

  zWidthMax_wvl = zAxEnd_wvl-zWTOriginOffset_wvl;
  zSVWidth_wvl = min(zSVWidthPre_wvl,  zWidthMax_wvl);
  
  xLen = ceil(xSVWidth_wvl/evParms.wt.xDelta_wvl);
  zLen = ceil(zSVWidth_wvl/evParms.wt.zDelta_wvl);

  PData(PDataIndWT).Origin = [xWTOriginOffset_wvl, 0, ...
                              zWTOriginOffset_wvl];

  evParms.wt.xSVWidth_wvl(ind) = xSVWidth_wvl;
  evParms.wt.zSVWidth_wvl(ind) = zSVWidth_wvl;

  evParms.wt.xEnd_wvl = xWTOriginOffset_wvl+...
      evParms.wt.xDelta_wvl*(xLen-1);   

  ax.xAx_wvl = ...
      xWTOriginOffset_wvl:evParms.wt.xDelta_wvl:evParms.wt.xEnd_wvl;

  evParms.wt.zEnd_wvl = zWTOriginOffset_wvl+...
      evParms.wt.zDelta_wvl*(zLen-1);   
  
  ax.zAx_wvl = ...
      zWTOriginOffset_wvl:evParms.wt.zDelta_wvl:evParms.wt.zEnd_wvl;  
  
end

% store these so we have them when we need to modify one at a time
wt.x = x;
wt.z = z;

PData(PDataIndWT).Size(1,1) = length(ax.zAx_wvl); % IQ samples for Doppler
PData(PDataIndWT).Size(1,2) = length(ax.xAx_wvl); % columns for Doppler
PData(PDataIndWT).Size(1,3) = 1; % single image page

PData(PDataIndWT).PDelta = [evParms.wt.xDelta_wvl, 0, ...
                            evParms.wt.zDelta_wvl];

ax.xAx_mm = ax.xAx_wvl*gate.lambda_mm;
ax.zAx_mm = ax.zAx_wvl*gate.lambda_mm;

[ax.X_mm, ax.Z_mm] = meshgrid(ax.xAx_mm, ax.zAx_mm);

evParms.wt.PDataIndWT = PDataIndWT;

evParms.wt.ax = ax;

%PData(PDataIndWT)

disp('leaving setpdatawtfn')
end

