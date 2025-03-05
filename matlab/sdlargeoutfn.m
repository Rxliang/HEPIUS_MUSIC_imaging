function cdiIm = sdlargeoutfn(varargin)

persistent cnt
persistent fifo
persistent fifoLen

persistent zeroIm

cdiIm = squeeze(varargin{1});

evParms = evalin('base', 'evParms');

if ~evParms.largeSDParms.SDLargeDisplayOverlay
  cdiIm=zeros(size(cdiIm));
  return
end

  
if isempty(cnt)
  cnt=0;
  fifoLen = 20; % frames
  fifo = zeros(1, fifoLen);
end

mx = maxall(cdiIm);

pos = mod(cnt, fifoLen)+1;

fifo(pos)=mx;

ind = find(fifo > eps);
  
%  fifoVals
if ~isempty(ind)
  medianMax = median(fifo(ind));
else
  medianMax=1;
end
  
cdiIm = cdiIm/medianMax*255;

cnt=cnt+1;
end

