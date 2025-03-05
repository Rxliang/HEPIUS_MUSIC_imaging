mfile = mfilename;
vtHome = getenv('MUSIC_HOME');
% mapping of adapter
inFile = [vtHome '/docs/pinmap_m5scd.mat'];
load(inFile, 'connMapAug', 'connMap');

ge408 = pinmap_ge408;

central168 = findcentraln(168,1:256); 

channels168 = ge408.channelNo(central168)

%return
central128 = findcentraln(150,1:256); 

selElements = 60+21:133+60-1;
% find the pin numbers on the 408 that the l8-18i will use

gePadsUsed = ge408.pad(selElements);

% now find the uta-260 connections for these
ind = [];
aperture =[];
for i = 1:length(gePadsUsed)
  ind = find(strcmp(connMap(1:128,2), gePadsUsed{i}));
  if ~isempty(ind)
    aperture(i) = ind;
    i
  end
  
end

aperture

GEPinNo=[];
for i = 1:length(central128)
  % this is the pin that this txd is connected to
  ind = find(ge408.channelNo == central128(i));
  GEPinNo(i)=ge408.channelNo(ind);
  
end

veraChannels = ge408.Connector1D(central128);
