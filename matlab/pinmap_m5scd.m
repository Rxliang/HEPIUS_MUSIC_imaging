mfile = mfilename;
ge408 = pinmap_ge408;

Trans.units = 'mm'; 
TransM5.name = 'GEM5ScD';
TransM5 = computeTrans_jon(TransM5);
 
% first 80 channels are at start of array
m5ElPerRow = TransM5.numelements/2;
m5CentralRowInd = 1:m5ElPerRow;
nChanAvailable = 64;
m5CentralRowCentral64Ind = findcentraln(nChanAvailable, m5CentralRowInd);

m5CentralRowOuterInd = setdiff(m5CentralRowInd, m5CentralRowCentral64Ind);

m5OuterRowInd = m5ElPerRow+1:2*m5ElPerRow;

outerElementsUsable = nChanAvailable*2-m5ElPerRow;
m5OuterRowCentralInd = findcentraln(outerElementsUsable, m5OuterRowInd);

ind128 = [m5CentralRowCentral64Ind m5CentralRowOuterInd ...
          m5OuterRowCentralInd];

% this vector now contains the  Verasonics logical indices of the
% channels we want, but we now need to remap to physical channels
% on the GE-408 connector

m5MapChannel128 = TransM5.Connector(ind128).';

m5PhysChannel128 = [];

% these should all be between 49 and 208

for i = 1:length(m5MapChannel128)
  m5PhysChannel128(i) = find(ge408.Connector1D==m5MapChannel128(i));
end

ind = find(m5PhysChannel128 > 208 | m5PhysChannel128 < 49);
if ~isempty(ind)
  error('Must use central channels on GE connector')
end

% check some channels

% E2 is 144 (144 is the transducer channel number)
% F2 is 173
%Trans.Connector = Trans.Connector([256 32]);
% Now 256 and 32 are the PHYSICAL CHANNELS EL_255 and EL_31, (E2
% and F2, resp)
 
if 0
  find(m5MapChannel128==144)
  find(m5MapChannel128==173)

  ge408.pad(256)
  ge408.pad(32)
end

% on M5Sc-D, our element 1 will be 165 (element 9) rather than 181
% (original)

%ge408.pad(m5PhysChannel128)
%ge408.channelNo(m5PhysChannel128)
%ge408.pad(TransM5.Connector)

%find(strcmp(ge408.pad(TransM5.Connector), 'E2'))
%find(strcmp(ge408.pad(TransM5.Connector), 'F2'))
%find(strcmp(ge408.pad(m5PhysChannel128), 'E2'))
%find(strcmp(ge408.pad(m5PhysChannel128), 'F2'))

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
connMap = [uta260sChannels.pinName' ge408.pad(m5PhysChannel128)' ...
           uta260sChannels.pinNumber'];

% mapping from UTA-260==verasonics physical, to ge channel numbers
connMapAug = [(1:128)' ge408.channelNo(m5PhysChannel128)'];

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


