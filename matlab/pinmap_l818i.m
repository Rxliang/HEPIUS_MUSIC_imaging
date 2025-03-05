% this is a 168-element transducer. We can use 128 elements
% maximum, in two apertures of 64 elements. These 128 elements
% can be chosen depending on mechanical design. Here, these are
% set to be the central 128 elements. The first aperture is the
% central 64, the second aperture is the outer 64

% the convention for GE connector is to use the central channels
% first, so GE physical channels 45-212 are connected to the 168
% elements of this probe

% we will make connections to GE physical channels 65-192
% (physChannel128InSeq below)

mfile = mfilename;
ge408 = pinmap_ge408;

TransL8.name = 'GEL8-18i';
Trans.name = 'GEL8-18i';
Trans.units = 'mm'; 
Trans = computeTrans_jon_4p0(Trans);

% transducer has 168 channels, we need to map 128 of these to the
% Cannon connector

numElements = length(Trans.Connector);
nChanAvailable = 128;
nChanPerMultiplex = 64;

centralInd128 = findcentraln(nChanAvailable, 1:numElements);
centralInd64 = findcentraln(nChanPerMultiplex, 1:numElements);

aperture1Ind = centralInd64;
aperture2Ind = setdiff(centralInd128, centralInd64)

ind128 = [aperture1Ind aperture2Ind]

% this vector now contains the  Verasonics logical indices of the
% channels we want, but we now need to remap to physical channels
% on the GE-408 connector

mapChannel128 = Trans.Connector(ind128).';

physChannel128InSeq = Trans.ConnectorES(centralInd128).';

physChannel128 = Trans.ConnectorES(ind128).';

% this is not needed, but is a check

L8PhysChannel128 = [];

% these should all be between 65 and 192

for i = 1:length(mapChannel128)
  L8PhysChannel128(i) = find(ge408.Connector1D==mapChannel128(i));
end

ind = find(L8PhysChannel128 > 192 | L8PhysChannel128 < 65);
if ~isempty(ind)
  error('Must use central channels on GE connector')
end

% check some channels

% ge408.pad{128} = 'K2' 
%ind = find(Trans.Connector==128)
%ge408.pad{4}
% E2 is 144 (144 is the transducer channel number)
% F2 is 173
%Trans.Connector = Trans.Connector([256 32]);
% Now 256 and 32 are the PHYSICAL CHANNELS EL_255 and EL_31, (E2
% and F2, resp)
 
if 1
%  find(mapChannel128==144)
  %find(mapChannel128==173)

  ge408.pad(256)
  ge408.pad(32)
end

uta260sOut = uta260s;

uta260sMap = [];

for i = 1:128
  searchStr = ['EL' num2str(i)];
  ind = find(strcmp(uta260sOut.pinName, searchStr));
  uta260sChannels.pinName{i} = searchStr;
  uta260sChannels.pinNumber{i} = uta260sOut.pinNumber{ind};
end

% signal name % GE % ATL

connMap = [];
connMap = [uta260sChannels.pinName' ge408.pad(physChannel128)' ...
           uta260sChannels.pinNumber'];

% mapping from UTA-260==verasonics physical, to ge channel numbers
connMapAug = [(1:128)' ge408.channelNo(physChannel128)'];

agnd.atl.ind = find(strcmp(uta260sOut.pinName, 'AGND'));

agnd.atl.pinNumber = uta260sOut.pinNumber(agnd.atl.ind);

agnd.atl.map = [repelem({'AGND'}, length(agnd.atl.ind),1) ...
                repelem({''}, length(agnd.atl.ind),1) ...
                agnd.atl.pinNumber'];

agnd.ge.mapPre = {'A1'; 'B1'; 'C1'; 'D1'; 'E1'; 'F1'; 'G1'; 'H1'; ...
               'J1'; 'K1'; 'L1'; 
               'A34'; 'C34'; 'D34'; 'E34'; 'F34'; 'G34'; 'H34'; ...
               'J34'; 'K34'; 'L34'; 'L22'; 'L23'; 'M16'; 'M17'};

spi.ge.sclk = 'L21';
ind = find(strcmp(uta260sOut.pinName, 'PCLK'));
spi.atl.pclk = uta260sOut.pinNumber{ind} %'A4'

spi.ge.sda = 'L20';
ind = find(strcmp(uta260sOut.pinName, 'PDAT'));
spi.atl.pdat = uta260sOut.pinNumber{ind};

ind = find(strcmp(uta260sOut.pinName, 'PPWR'));
spi.atl.ppwr = uta260sOut.pinNumber{ind};

ind = find(strcmp(uta260sOut.pinName, 'PWEN'));
spi.atl.pwen = uta260sOut.pinNumber{ind};

%led.ge.ppg = 'M2';
pwr.ge.vcc = {'M23' 'M24'};
pwr.ge.v12vp = 'M18';

mux.ge.hvp = {'M33', 'M34'};
mux.ge.hvn = {'M29', 'M30'};


ind = find(strcmp(uta260sOut.pinName, 'VML'));
mux.atl.vml = uta260sOut.pinNumber{ind};

ind = find(strcmp(uta260sOut.pinName, 'VPP'));
mux.atl.vpp = uta260sOut.pinNumber{ind};

ind = find(strcmp(uta260sOut.pinName, 'VNN'));
mux.atl.vnn = uta260sOut.pinNumber{ind};

% route PPWR to 5V on GE

% VML This is the logic supply for the HV Mux chips. This is dependent ...
%      ion the type of chip used and vary from 3.3V to ~12V.
 
% use VML to source 12V and 5V.

% bank B
banks = 'BD';
for j = 1:length(banks)
  for i = 2:2:32
    agnd.ge.mapPre = [agnd.ge.mapPre; [banks(j) num2str(i)]];
  end
end

banks = 'GJ';
for j = 1:length(banks)
  for i = 3:2:33
    agnd.ge.mapPre = [agnd.ge.mapPre; [banks(j) num2str(i)]];
  end
end

agnd.ge.map = [repelem({'AGND'}, length(agnd.ge.mapPre),1) ...
               agnd.ge.mapPre ...
               repelem({''}, length(agnd.ge.mapPre),1)];

connMap = [connMap; agnd.ge.map; agnd.atl.map];

connMap = [connMap; 
           {'SCLK', spi.ge.sclk, spi.atl.pclk};
           {'SDAT', spi.ge.sda, spi.atl.pdat};
           {'VCC', pwr.ge.vcc{1}, spi.atl.ppwr};
           {'VCC', pwr.ge.vcc{2}, spi.atl.ppwr}; 
           %{'PPG', led.ge.ppg, spi.atl.pwen};
           {'HVP', mux.ge.hvp{1}, mux.atl.vpp};
           {'HVP', mux.ge.hvp{2}, mux.atl.vpp};
           {'HVN', mux.ge.hvn{1}, mux.atl.vnn};
           {'HVN', mux.ge.hvn{2}, mux.atl.vnn};           
           {'V12VP', pwr.ge.v12vp, mux.atl.vml}];

vtHome = getenv('MUSIC_HOME');
outFile = [vtHome '/docs/' mfile '.xlsx'];
outMat = [vtHome '/docs/' mfile '.mat'];
unix(['rm ' outFile]);

connTab = cell2table(connMap);
connTab.Properties.VariableNames = {'SignalName', 'DLP408', 'DL260'};

writetable(connTab, outFile);
lslrt(outFile);

save(outMat, 'connTab', 'connMap', 'connMapAug');
lslrt(outMat);


