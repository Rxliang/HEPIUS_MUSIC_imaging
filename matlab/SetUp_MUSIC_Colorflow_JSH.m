%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Last updated 2022-May-18 by Jonah (joharmon@uw.edu)
%
% This script defaults to live wide beam Bmode imaging for alignment and
% identification of an imaging plane of interest, with a GUI control to
% fire a singleshot plane wave Doppler acquisition. Data will be
% reconstructed and displayed immediately following acquisition -
% specifically, flow power and flow velocity images. This has been tested
% on Ubuntu 18.04, kernel 4.18.0-25-generic.
%
% A number of controls have been built into the script to provide
% flexibility for the end user with regards to data acquisition,
% processing, and saving (e.g., raw IQ, processed IQ, or end results).
% These are aggregated near the top of the script for easy access.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all

% Define your Vantage directory path
%vantage_directory = '~/Documents/Vantage-4.6.2-2110271004';
%vantage_directory = '~/Documents/Vantage-4.6.2-2110271004';
vantage_directory = 'c:/Users/223022005/Vantage-4.9.2-2308102000/';

%% Setup
cd(vantage_directory)
activate

%% User-defined options
% What data to save and where to save it
opts.save_path = '~/Documents/test_save/';  % Path for saving data - ideally on SSD for speed
opts.save_RF = 1;                           % Raw RF data - SLOW if set to 1
opts.save_IQ = 1;                           % Unprocessed IQ ensemble - SLOW if set to 1
opts.save_Doppler_output = 1;               % .mat file with raw flow power and velocity images
opts.save_images = 1;                       % Flow power and velocity images - .fig and .png
opts.save_bmode_images = 1;                 % Single Bmode - .png and .fig

% Doppler data acquisition, processing, and display
opts.singleshot_PRF = 2000;                           % Effective PRF, in Hz - will dictate Nyquist velocity
opts.singleshot_PRI = 600;                            % Number of frames in Doppler ensemble
opts.num_angles = 5;                                  % Number of angles for coherent compounding
opts.filter_type = 'SVD';                             % 'SVD' only for now
opts.filter_thresholds = [50 opts.singleshot_PRI];    % For SVD wall filtering
opts.power_display = [135 180];                        % In dB
opts.velocity_display = [-3 3];                       % In cm/sec
opts.power_processing = 'R1';                         % 'Cleaner' or 'R1'
opts.velocity_power_threshold = 120;                  % For thresholding vessels vs. background

% Live Bmode
opts.num_cine_frames = 50;                            % Number of DisplayWindow frames for cine loop

% Bmode display alongside CF singleshot
opts.bmode_display_dynamic_range = 60;                % Display dynamic range in dB 

%% Define Bmode parameters for live imaging
P.numTx = 32;            % no. of elements in TX aperture.
P.numRays = 48;          % no. of rays in frame
P.txFocusMm = 20;        % transmit focal pt 
P.dtheta = 10*(pi/180);  % angle delta for beams steered left or right
P.startDepthMm = 0.0;    % startDepth in mm
P.endDepthMm = 12.0;     % endDepth in mm
P.maxDepthMm = 20.0;     % maxDepth for RangeChange and RcvBuffer

%% Define colorflow parameters for singleshot
% Using CF struct to keep independent of Bmode parameters in P struct
CF.startDepthMm = 0.0;           % startDepth in mm
CF.endDepthMm = 12.0;            % endDepth in mm
CF.maxDepthMm = 20.0;            % maxDepth for RangeChange and RcvBuffer
CF.opening_angle = (24*pi/180);  % 24 degree sweep for angle compounding

% Defining transmit angles
na = opts.num_angles;         % Number of angles, coherently compounded
if (na > 1)
    CF.dtheta = CF.opening_angle/(na-1); 
    CF.startAngle = -CF.opening_angle/2; 
else
    CF.dtheta = 0; CF.startAngle=0;
end 

%% Define general parameters
TWFreq = 15.625;              % Center frequency
sample_Mode_b = 'BS100BW';    % 100% bandwidth
samplesperwave_b = 2; 
TimeTagEna = 0;               % Initial state of time tag function
interp_factor = 1;            % No interpolation

%% Specify Resource.Parameters
Resource.Parameters.numTransmit = 64;     % Number of transmit channels.
Resource.Parameters.numRcvChannels = 64;  % Number of receive channels.
Resource.Parameters.speedOfSound = 1540;  % Speed of sound in m/s
Resource.Parameters.verbose = 2;          % For debug
Resource.Parameters.initializeOnly = 0;   
Resource.Parameters.simulateMode = 0;     % Set to 1 to run in simulation

%% Specify Transducer
Trans.name = 'MUSIC';            % 'MUSIC'; 
Trans.units = 'wavelengths';
Trans = computeTrans2(Trans);
Trans.maxHighVoltage = 25;       % Set maximum high voltage limit
Trans.id = 305746;
Trans.name = 'custom';
Trans.frequency = 125./9;

% Determine required buffer size
demodFreq = TWFreq;                        % Demodulation frequency
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
P.txFocus = P.txFocusMm*scaleToWvl;
P.startDepth = P.startDepthMm*scaleToWvl;  % startDepth in wavelength
P.endDepth = P.endDepthMm*scaleToWvl;
P.maxDepth = P.maxDepthMm*scaleToWvl;
maxBufLength = ceil(sqrt(P.maxDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
maxBufSizePerAcq = 64*ceil(maxBufLength*8*(demodFreq/Trans.frequency)/64);

% Compute maximal memory required at max depth -- colorflow
maxDepth_wl = CF.maxDepthMm*Trans.frequency/1540*1000;
maxAcqLength = ceil(sqrt(maxDepth_wl^2 + ((Trans.numelements-1)*Trans.spacing)^2));
wlsPer128_b = 128/(samplesperwave_b); % wavelengths in 128 samples for 4 samplesPerWave
required_memory_b = na*2*wlsPer128_b*(ceil(maxAcqLength/wlsPer128_b))*samplesperwave_b;

%% Specify Bmode PData
% Recon pixel grid information, Bmode
PData(1).PDelta = [Trans.spacing/4, 0, 1/4];
opts.mmPerPxZ = PData(1).PDelta(3) / scaleToWvl;
opts.mmPerPxX = PData(1).PDelta(1) / scaleToWvl;
PData(1).Size(1,1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3));
PData(1).Size(1,2) = ceil((Trans.numelements*Trans.spacing)/PData(1).PDelta(1));
PData(1).Size(1,3) = 1;
PData(1).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of uppr lft crnr.

% Predefine all the Region structures, then modify them below.
PData(1).Region = repmat(struct('Shape',struct('Name','Parallelogram',...
                                               'Position',[0,0,P.startDepth],...
                                               'width',5*Trans.spacing*Trans.numelements/(P.numRays),...
                                               'height',P.endDepth-P.startDepth,...
                                               'angle', 0.0)),1,3*P.numRays); % default is no steering
                                           
% Compute the x coords of the TX beam centers
TxOrgX = (-31.5*Trans.spacing):(63*Trans.spacing/(P.numRays-1)):(31.5*Trans.spacing);

% Specify P.numRays rectangular regions centered on TX beam origins (use default angle of 0.0).
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
m = P.numRays;

% Define numRays steered left parallelogram regions, centered on TX beam 
% origins. Adjust the angle so that the steering goes to zero over 8 beams 
% at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = -((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = -((P.numRays-n)/8)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
m = m + P.numRays;

% Define numRays steered right parallelogram regions, centered on TX beam 
% origins. Adjust the angle so that the steering goes to zero over 8 beams 'Units', 'normalized', 'Position', [0.1, 0.5, 0.25, 0.25], ...
% at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = ((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = ((P.numRays-n)/8)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
PData(1).Region = computeRegions(PData(1));

%% Specify colorflow singleshot PData
% No PData.Region specified, so a default Region for the entire PData array 
% will be created by computeRegions.
% Same pixels dimensions as Bmode; same imaging depth as Bmode
PData(2).PDelta = [Trans.spacing/4, 0, 1/4];
opts.mmPerPxZ_dop = PData(2).PDelta(3) / scaleToWvl;
opts.mmPerPxX_dop = PData(2).PDelta(1) / scaleToWvl;
PData(2).Size(1) = ceil((P.endDepth-P.startDepth)/PData(2).PDelta(3));
PData(2).Size(2) = ceil((Trans.numelements*Trans.spacing)/PData(2).PDelta(1));
PData(2).Size(3) = 1;
PData(2).Origin = [-Trans.spacing*(Trans.numelements-1)/2,0,P.startDepth]; % x,y,z of upper lft crnr.

%% Specify Media object
pt1;
Media.attenuation = -0.5;
Media.function = 'movePoints';

%% Specify Resources
% Rows per frame used here allows for all rays, with range of up to 400 
% wavelengths.
% Receive buffers
num_rcv_frames = 12;
Resource.RcvBuffer(1).datatype = 'int16';
Resource.RcvBuffer(1).rowsPerFrame = P.numRays*maxBufSizePerAcq; 
Resource.RcvBuffer(1).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(1).numFrames = num_rcv_frames;
Resource.RcvBuffer(2).datatype = 'int16';
Resource.RcvBuffer(2).rowsPerFrame = P.numRays*maxBufSizePerAcq;
Resource.RcvBuffer(2).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(2).numFrames = num_rcv_frames; 
Resource.RcvBuffer(3).datatype = 'int16';
Resource.RcvBuffer(3).rowsPerFrame = P.numRays*maxBufSizePerAcq; 
Resource.RcvBuffer(3).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(3).numFrames = num_rcv_frames;

% Inter buffer
Resource.InterBuffer(1).datatype = 'complex';
Resource.InterBuffer(1).numFrames = 1;  % One intermediate buffer needed.

% Image buffer
Resource.ImageBuffer(1).datatype = 'double';
Resource.ImageBuffer(1).numFrames = 10;

%% Colorflow singleshot buffers
% Receive
Resource.RcvBuffer(4).datatype = 'int16';
Resource.RcvBuffer(4).rowsPerFrame = required_memory_b*opts.singleshot_PRI;
Resource.RcvBuffer(4).colsPerFrame = Resource.Parameters.numRcvChannels;
Resource.RcvBuffer(4).numFrames = 1;  

% Inter buffer (for long, single acquisition Doppler reconstruction)
Resource.InterBuffer(2).datatype = 'complex';
Resource.InterBuffer(2).numFrames = 1;
Resource.InterBuffer(2).pagesPerFrame = opts.singleshot_PRI;

%% Display - Bmode only, custom display for colorflow singleshot
Resource.DisplayWindow(1).Title = 'MUSIC Bmode Live Imaging';
Resource.DisplayWindow(1).pdelta = 0.25;
ScrnSize = get(0,'ScreenSize');
DwWidth = ceil(PData(1).Size(2)*PData(1).PDelta(1)/Resource.DisplayWindow(1).pdelta);
DwHeight = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);
Resource.DisplayWindow(1).Position = [250,(ScrnSize(4)-(DwHeight+150))/2, ...  % lower left corner position
                                      DwWidth, DwHeight];
Resource.DisplayWindow(1).ReferencePt = [PData(1).Origin(1),0,PData(1).Origin(3)];   % 2D imaging is in the X,Z plane
Resource.DisplayWindow(1).Type = 'Matlab'; % Can specify 'Verasonics' if desired
Resource.DisplayWindow(1).numFrames = opts.num_cine_frames;
Resource.DisplayWindow(1).AxesUnits = 'mm';
Resource.DisplayWindow(1).Colormap = gray(256);

%% TPC
% Monitoring voltage, Bmode
TPC(1).name = 'Monitoring'; 
TPC(1).maxHighVoltage = 25;

% Doppler acquisition voltage
TPC(2).name = 'Doppler';    
TPC(2).maxHighVoltage = 25;

%% Specify Transmit waveform structure
% Bmode
TW(1).type = 'parametric';
TW(1).Parameters = [Trans.frequency,0.97,1,1]; 

% Doppler
TW(2).type = 'parametric';
TW(2).Parameters = [15.625,.67,4,1]; % Two cycle pulse

%% Specify TX structure array
% Bmode plus na transmits for colorflow
TX = repmat(struct('waveform', 1, ...
                   'Origin', [0.0,0.0,0.0], ...
                   'Apod', zeros(1,Resource.Parameters.numTransmit), ...
                   'focus', P.txFocus, ...
                   'Steer', [0.0,0.0], ...
                   'Delay', zeros(1,Resource.Parameters.numTransmit), ...
                   'TXPD', [], ...
                   'peakCutOff', 0.5, ...
                   'peakBLMax', 15.0), 1, 3*P.numRays + na);

if strcmp(Trans.units, 'wavelengths')
    scaleToWvl = 1;
end

%% Set event specific TX attributes
for n = 1:P.numRays
    TX(n).Origin(1) = TxOrgX(n);
    
    % Compute transmit aperture apodization
    TX(n).Apod = +(((scaleToWvl*Trans.ElementPos(:,1))>(TxOrgX(n)-Trans.spacing*P.numTx/2))& ...
                 ((scaleToWvl*Trans.ElementPos(:,1))<(TxOrgX(n)+Trans.spacing*P.numTx/2)))';
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
end
m = P.numRays;

for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [-((n-1)/8)*P.dtheta,0.0];
    elseif n>(P.numRays-8)
        TX(n+m).Steer = [-((P.numRays-n)/8)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [-P.dtheta,0.0];
    end
    
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end
m = m + P.numRays;

for n = 1:P.numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [((n-1)/8)*P.dtheta,0.0];
    elseif n>(P.numRays-8)
        TX(n+m).Steer = [((P.numRays-n)/8)*P.dtheta,0.0];
    else
        TX(n+m).Steer = [P.dtheta,0.0];
    end
    
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
end
start_CF_TX = n+m+1;

% % Calculate TXPD
% fprintf('Computing Bmode TXPD...\n')
% steps = 3*P.numRays;
% for i = 1:steps
%     TX(i).TXPD = computeTXPD(TX(i),PData(1));
% end

%% Colorflow transmits
% Set defaults first
for i = start_CF_TX:1:size(TX, 2)
    TX(i).waveform = 2;
    TX(i).Origin = [0.0, 0.0, 0.0];
    TX(i).focus = 0.0;
    TX(i).Steer = [0.0, 0.0];
    TX(i).Delay = zeros(1,Resource.Parameters.numTransmit);
    TX(i).Apod = kaiser(Resource.Parameters.numTransmit,1)';
    TX(i).peakCutOff = 0.5; % Need to change this? From Bmode
    TX(i).peakBLMax = 15.0; % Need to change this? From Bmode
end

% Define actual steered plane wave transmits
counter = start_CF_TX-1;
for n = 1:na
    angle = CF.startAngle+(n-1)*CF.dtheta;
    TX(n+counter).Steer = [angle,0.0];
    TX(n+counter).Apod= kaiser(Resource.Parameters.numTransmit,1)';
    TX(n+counter).Delay = computeTXDelays(TX(n+counter));
end

%% Specify Receive structure arrays
% Bmode + na*PRI for colorflow
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
Receive = repmat(struct('Apod', ones(1,Trans.numelements), ...
                        'startDepth', P.startDepth, ...
                        'endDepth', maxAcqLength, ...
                        'TGC', 1, ...
                        'bufnum', 1, ...
                        'framenum', 1, ...
                        'acqNum', 1, ...
                        'sampleMode', 'NS200BW', ...
                        'mode', 0, ...
                        'callMediaFunc', 0), 1, 3*P.numRays*Resource.RcvBuffer(1).numFrames + na*opts.singleshot_PRI);
                    
RcvProfile(1).LnaGain = 18;   % Profile used for Bmode live imaging
RcvProfile(2).LnaGain = 24;   % Profile used for Doppler
RcvProfile(2).LnaZinSel = 31; % Force high-Z state for best Doppler sensitivity  

%% Define filters - MAY NEED TO ADJUST FOR MUSIC
noLPF = [zeros(1,11) 1];

wideinput = [+0.00000 +0.00000 +0.00131 +0.00000 +0.00525 +0.00000 +0.01172 ...
             +0.00000 +0.01620 +0.00000 +0.00891 +0.00000 -0.02042 +0.00000 ...
             -0.07434 +0.00000 -0.14124 +0.00000 -0.19757 +0.00000 +0.78033];

%% Set event specific Receive attributes for each frame
% Bmode
Receive(1).callMediaFunc = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    k = 3*P.numRays*(i-1);
    Receive(k+1).callMediaFunc = 1;
    for j = 1:P.numRays
        Receive(k+j).framenum = i;
        Receive(k+j).acqNum = j;
        Receive(k+P.numRays+j).bufnum = 2;
        Receive(k+P.numRays+j).framenum = i;
        Receive(k+P.numRays+j).acqNum = j;
        Receive(k+2*P.numRays+j).bufnum = 3;
        Receive(k+2*P.numRays+j).framenum = i;
        Receive(k+2*P.numRays+j).acqNum = j;
    end
end
LastRxNum = k+2*P.numRays+j;

% Colorflow singleshot
k = LastRxNum;
for j = 1:na*opts.singleshot_PRI % Loop over ensemble length
    Receive(k+j).sampleMode = sample_Mode_b;
    Receive(k+j).Apod(1:Trans.numelements) = 1;
    Receive(k+j).framenum = 1;
    Receive(k+j).acqNum = j;
    Receive(k+j).TGC = 2;
    Receive(k+j).bufnum = 4; % Need to refer to proper buffer
    Receive(k+j).LowPassCoef = noLPF;
    Receive(k+j).InputFilter = wideinput;
end

%% Specify TGC Waveform structures
% Bmode TGC
TGC(1).CntrlPts = [0,225,410,475,549,659,736,791];
TGC(1).rangeMax = P.endDepth;
TGC(1).Waveform = computeTGCWaveform(TGC(1));

% Doppler TGC
TGC(2).CntrlPts = [690,1023,1023,1023,1023,1023,1023,1023];
TGC(2).rangeMax = P.endDepth;
TGC(2).Waveform = computeTGCWaveform(TGC(2));

%% Specify Recon structure arrays
% We need three Recon structures, one for each steering direction, then one
% for colorflow singleshot
Recon = repmat(struct('senscutoff', 0.5, ...
               'pdatanum', 1, ...
               'IntBufDest', [1,1], ...
               'ImgBufDest', [1,-1], ...
               'rcvBufFrame',-1, ...
               'RINums',1:P.numRays), 1, 3+1);

% Set specific Recon attributes.
Recon(2).RINums = (P.numRays+1):(2*P.numRays);
Recon(3).RINums = (2*P.numRays+1):(3*P.numRays);
end_bmode_RI = 3*P.numRays;

% Colorflow singleshot
Recon(4).RINums = end_bmode_RI+1:end_bmode_RI + na*opts.singleshot_PRI;
Recon(4).IntBufDest = [2,1]; 
Recon(4).ImgBufDest = [0,0];
Recon(4).pdatanum = 2;
total_RI_nums = end_bmode_RI + na*opts.singleshot_PRI;

%% Define ReconInfo structures
ReconInfo = repmat(struct('mode', 'accumIQ', ...  % default is to accumulate IQ data.
                   'Pre',[],...
                   'Post',[],...
                   'txnum', 1, ...
                   'rcvnum', 1, ...
                   'pagenum', 1, ...
                   'scaleFactor', 0.5, ...
                   'regionnum', 0), 1, total_RI_nums);
               
% Change defaults for colorflow RIs
for i = end_bmode_RI+1:end_bmode_RI + na*opts.singleshot_PRI
    ReconInfo(i).scaleFactor = 0.5; % Change this to 0.5 to match Bmode?
    ReconInfo(i).regionnum = 1;
end

%% Set specific ReconInfo attributes
ReconInfo(1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j).txnum = j;
    ReconInfo(j).rcvnum = j;
    ReconInfo(j).regionnum = j;
end
ReconInfo(P.numRays).Post = 'IQ2IntensityImageBuf';
k = P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = j+k;
    ReconInfo(j+k).rcvnum = j+k;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';
k = k + P.numRays;
ReconInfo(k+1).Pre = 'clearInterBuf';
for j = 1:P.numRays
    ReconInfo(j+k).txnum = j+k;
    ReconInfo(j+k).rcvnum = j+k;
    ReconInfo(j+k).regionnum = j+k;
end
ReconInfo(k+P.numRays).Post = 'IQ2IntensityImageBuf';

% Colorflow singleshot - loop through angles within each page
for k=1:opts.singleshot_PRI
    for j = 1:na
        ReconInfo(end_bmode_RI + (k-1)*na + j).txnum = start_CF_TX-1 + j;
        ReconInfo(end_bmode_RI + (k-1)*na + j).rcvnum = LastRxNum + (k-1)*na + j;
        ReconInfo(end_bmode_RI + (k-1)*na + j).pagenum = k;
        ReconInfo(end_bmode_RI + (k-1)*na + j).regionnum = 1;
        if j == 1
            ReconInfo(end_bmode_RI + (k-1)*na + j).mode = 'replaceIQ';
        elseif j == na
            ReconInfo(end_bmode_RI + (k-1)*na + j).mode = 'accumIQ';
        else
            ReconInfo(end_bmode_RI + (k-1)*na + j).mode = 'accumIQ';
        end
    end
end

%% Specify Process structure array
% Image display for Bmode live
pers = 20;
Process(1).classname = 'Image';
Process(1).method = 'imageDisplay';
Process(1).Parameters = {'imgbufnum',1,...   % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum',1,...    % number of PData structure to use
                         'pgain',3.0,...     % pgain is image processing gain
                         'reject',2,...      % reject level
                         'persistMethod','simple',...
                         'persistLevel',pers,...
                         'interpMethod','4pt',...
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','runAverage3',...
                         'compressMethod','log',...
                         'compressFactor',50,...
                         'mappingMethod','full',...
                         'display',1,...      % display image after processing
                         'displayWindow',1};
                   
% Save raw RF data if flag is set
Process(2).classname = 'External';
Process(2).method = 'save_raw_RF';
Process(2).Parameters = {'srcbuffer', 'receive', ...
                         'srcbufnum', 4, ...
                         'dstbuffer', 'none'};

% Save IQ ensemble if flag is set
Process(3).classname = 'External';
Process(3).method = 'save_IQ_ensemble';
Process(3).Parameters = {'srcbuffer', 'inter', ...
                         'srcbufnum', 2, ...
                         'dstbuffer', 'none'};

% Doppler processing post-recon - display and save rolled into this
Process(4).classname = 'External';
Process(4).method = 'process_Doppler';
Process(4).Parameters = {'srcbuffer', 'inter', ...
                         'srcbufnum', 2, ...
                         'dstbuffer', 'none'};
                     
% Bmode frame display and save - grab one to match each CF acq
Process(5).classname = 'External';
Process(5).method = 'save_Bmode';
Process(5).Parameters = {'srcbuffer', 'image', ...
                         'srcbufnum', 1, ...
                         'dstbuffer', 'none'};
                     
% Reset start event                
Process(6).classname = 'External';
Process(6).method = 'reset_start_event';
Process(6).Parameters = {'srcbuffer', 'none', 'dstbuffer', 'none'};

%% Specify SeqControl structure arrays
nsc = 1;

% Jump to first event
sc_jump = nsc;
SeqControl(nsc).command = 'jump';
SeqControl(nsc).argument = 1;
nsc = nsc+1;

% Time between synthetic aperture acquisitions
sc_ttna1 = nsc;
SeqControl(nsc).command = 'timeToNextAcq'; 
SeqControl(nsc).argument = 220;  % usec
nsc = nsc+1;

% Optional time between frames (not needed for synchronous operation)
sc_ttna2 = nsc;
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = 30000;  % 30000 usec = 30msec time between frames
nsc = nsc+1;

% Return to Matlab
sc_return = nsc;
SeqControl(nsc).command = 'returnToMatlab';
nsc = nsc+1;

%% Colorflow-specific SeqControl
% Change to Profile 2 (Doppler)
sc_dop_tpc = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).condition = 'immediate';
SeqControl(nsc).argument = 2;
nsc = nsc+1;

% Change to Profile 1 (Bmode)
sc_bmode_tpc = nsc;
SeqControl(nsc).command = 'setTPCProfile';
SeqControl(nsc).condition = 'immediate';
SeqControl(nsc).argument = 1;
nsc = nsc+1;

% PRF for plane wave colorflow - equal inter-frame and inter-angle PRF
sc_dop_ttna = nsc;
SeqControl(nsc).command = 'timeToNextAcq';
SeqControl(nsc).argument = round(1/(na*opts.singleshot_PRF*1e-06));
nsc = nsc+1;

% Set receive profile 1
sc_rcv_profile_1 = nsc;
SeqControl(nsc).command = 'setRcvProfile';
SeqControl(nsc).argument = 1;
nsc = nsc+1;

% Set receive profile 2
sc_rcv_profile_2 = nsc;
SeqControl(nsc).command = 'setRcvProfile';
SeqControl(nsc).argument = 2;
nsc = nsc+1;

% Info regarding max sampled velocity at specified PRF
disp([num2str(opts.singleshot_PRI*((na)*SeqControl(sc_dop_ttna).argument)/1000),...
    ' ms singleshot acquisition time; ', num2str(opts.singleshot_PRF), ' Hz effective PRF; ',...
    num2str(100*(Resource.Parameters.speedOfSound*opts.singleshot_PRF)/(4*TWFreq*10^6)),...
    ' cm/sec max sampled velocity'])
opts.max_velocity = 100*(Resource.Parameters.speedOfSound*opts.singleshot_PRF)/(4*TWFreq*10^6); % cm/sec

%% Specify Event sequence - live Bmode
n = 1;
for i = 1:Resource.RcvBuffer(1).numFrames
    % Acquire frame with unsteered wide beams
    k = 3*P.numRays*(i-1);
    for j = 1:P.numRays       
        Event(n).info = 'Acquire aperture';
        Event(n).tx = j;
        Event(n).rcv = k+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = sc_ttna1;
        
        % Set rcv profile if on first transmit
        if i*j == 1
            Event(n).seqControl = [sc_rcv_profile_2, sc_ttna1];
        end
        
        n = n+1;
    end
    
    % Modify last Event's SeqControl for ttna, transfer frame to host
    Event(n-1).seqControl = [sc_ttna2, nsc];
      SeqControl(nsc).command = 'transferToHost';
      nsc = nsc+1;
    
    Event(n).info = 'Recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;Process(4).Parameters = {'srcbuffer', 'inter', ...
                         'srcbufnum', 2, ...
                         'dstbuffer', 'none'};
    Event(n).recon = 1;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered left wide beams
    for j = 1:P.numRays      
        Event(n).info = 'Acquire aperture';
        Event(n).tx = j+P.numRays;
        Event(n).rcv = k+P.numRays+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = sc_ttna1;
        n = n+1;
    end
    
    % Modify last Event's SeqControl for ttna, transfer frame to host
    Event(n-1).seqControl = [sc_ttna2, nsc]; 
      SeqControl(nsc).command = 'transferToHost'; 
      nsc = nsc+1;

    Event(n).info = 'Recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 2;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    n = n+1;

    % Acquire frame with steered right wide beams
    for j = 1:P.numRays       
        Event(n).info = 'Acquire aperture';
        Event(n).tx = j+2*P.numRays;
        Event(n).rcv = k+2*P.numRays+j;
        Event(n).recon = 0;
        Event(n).process = 0;
        Event(n).seqControl = sc_ttna1;
        n = n+1;
    end
    
    % Modify last Event's SeqControl for ttna, transfer frame to host
    Event(n-1).seqControl = [sc_ttna2, nsc];
      SeqControl(nsc).command = 'transferToHost'; 
      nsc = nsc+1;

    Event(n).info = 'Recon and process';
    Event(n).tx = 0;
    Event(n).rcv = 0;
    Event(n).recon = 3;
    Event(n).process = 1;
    Event(n).seqControl = 0;
    if i ~= Resource.RcvBuffer(1).numFrames  % Exit to Matlab every 3rd frame
        Event(n).seqControl = sc_return;
    end
    n = n+1;
end

% Loop to continue live Bmode imaging
Event(n).info = 'Jump to Event 1';
Event(n).tx = 0;
Event(n).rcv = 0;
Event(n).recon = 0;
Event(n).process = 0;
Event(n).seqControl = sc_jump;
n = n+1;

%% Colorflow singleshot
singleshot_start = n;

% Switch transmit power controller
Event(n).info = 'Switch to TPC2 and wait';
Event(n).tx = 0;         % no tx
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no recon
Event(n).process = 0;    % no processing
Event(n).seqControl = sc_dop_tpc;
n = n+1;

% Acquire all frames defined in RcvBuffer
for i = 1:opts.singleshot_PRI
    k = na*(i-1)+LastRxNum;
    for j = 1:na
        Event(n).info = 'Acquire plane wave for Colorflow';
        Event(n).tx = start_CF_TX-1 + j;   % use plane wave tx
        Event(n).rcv = k+j;                % use plane wave rcvs
        Event(n).recon = 0;                % no reconstruction.
        Event(n).process = 0;              % no processing
        Event(n).seqControl = sc_dop_ttna; % time between syn. aper. acqs.
        if i*j == 1
            Event(n).seqControl = [sc_rcv_profile_2, sc_dop_ttna];
        end
        n = n+1;
    end
end

% Only transfer buffer at end of Doppler sequence to optimize data transfer
Event(n-1).seqControl = nsc; % use SeqControl structs defined below.
SeqControl(nsc).command = 'transferToHost';
lastTTHnsc=nsc;
nsc = nsc + 1;

% Make sure last transfer has cleared before trying to reconstruct
Event(n).info = 'waitForTransferComplete';
Event(n).tx = 0;            % no tx
Event(n).rcv = 0;           % no rcv
Event(n).recon = 0;         % no recon
Event(n).process = 0;       % no process
Event(n).seqControl = nsc;
SeqControl(nsc).command = 'waitForTransferComplete';
SeqControl(nsc).argument = lastTTHnsc;
nsc = nsc + 1;
n = n+1;

Event(n).info = 'Reconstruct';
Event(n).tx = 0;            % no TX
Event(n).rcv = 0;           % no Rcv
Event(n).recon = 4;         % singleshot recon
Event(n).process = 0;       % no process structure
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Save RF if flagged';
Event(n).tx = 0;            % no TX
Event(n).rcv = 0;           % no Rcv
Event(n).recon = 0;
Event(n).process = 2;       
Event(n).seqControl = 0;
n = n+1;

Event(n).info = 'Save IQ if flagged';
Event(n).tx = 0;                % no tx
Event(n).rcv = 0;               % no rcv
Event(n).recon = 0;             % no recon
Event(n).process = 3;         
Event(n).seqControl = 0;        % no instruction
n = n+1;

Event(n).info = 'Doppler process, display, save';
Event(n).tx = 0;                % no tx
Event(n).rcv = 0;               % no rcv
Event(n).recon = 0;             % no recon
Event(n).process = 4;           
Event(n).seqControl = 0;        % no instruction
n = n+1;

Event(n).info = 'Bmode display, save';
Event(n).tx = 0;                % no tx
Event(n).rcv = 0;               % no rcv
Event(n).recon = 0;             % no recon
Event(n).process = 5;           
Event(n).seqControl = 0;        % no instruction
n = n+1;

% Reset start event to live imaging, don't want a new burst acq every time
% we unfreeze
Event(n).info = 'Reset start event to monitoring loop';
Event(n).tx = 0;         % no tx
Event(n).rcv = 0;        % no rcv
Event(n).recon = 0;      % no recons
Event(n).process = 6;    
Event(n).seqControl = 0; % no instruction
n = n+1;

% Switch transmit power controller back to Bmode
Event(n).info = 'Switch to TPC1';
Event(n).tx = 0;                     % no tx
Event(n).rcv = 0;                    % no rcv
Event(n).recon = 0;                  % no recon
Event(n).process = 0;                % no processing
Event(n).seqControl = sc_bmode_tpc;  % TPC1
n = n+1;

% Quit to matlab
Event(n).info = 'Return to Matlab';
Event(n).tx = 0;                        % no tx
Event(n).rcv = 0;                       % no rcv
Event(n).recon = 0;                     % no recon
Event(n).process = 0;                 
Event(n).seqControl = sc_return;        % Stop executing instructions
n = n+1;

%% User specified UI Control Elements
% Sensitivity Cutoff
UI(1).Control =  {'UserB7','Style','VsSlider','Label','Sens. Cutoff',...
                  'SliderMinMaxVal',[0,1.0,Recon(1).senscutoff],...
                  'SliderStep',[0.025,0.1],'ValueFormat','%1.3f'};
UI(1).Callback = text2cell('%SensCutOffCallback');

% Range Change
scaleToWvl = 1;
AxesUnit = 'mm';
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if ~strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        AxesUnit = 'wls';
        scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);
    end
end
UI(2).Control = {'UserA1','Style','VsSlider','Label',['Range (',AxesUnit,')'],...
                 'SliderMinMaxVal',[2,P.maxDepthMm,P.endDepthMm]*scaleToWvl,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(2).Callback = text2cell('%RangeChangeCallback');

% Transmit focus change
UI(3).Control = {'UserB4','Style','VsSlider','Label',['TX Focus (',AxesUnit,')'],...
                 'SliderMinMaxVal',[10,600,P.txFocusMm]*scaleToWvl,'SliderStep',[0.1,0.2],'ValueFormat','%3.0f'};
UI(3).Callback = text2cell('%TxFocusCallback');

% TX Aperture change
UI(4).Control = {'UserB5','Style','VsSlider','Label','TX Aper',...
                 'SliderMinMaxVal',[1,128,P.numTx],'SliderStep',[1/128,1/32],'ValueFormat','%3.0f'};
UI(4).Callback = text2cell('%ApertureCallback');

% Peak CutOff
UI(5).Control = {'UserB2','Style','VsSlider','Label','Peak Cutoff',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakCutOff],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(5).Callback = text2cell('%PeakCutOffCallback');

% Max. Burst Length
UI(6).Control = {'UserB1','Style','VsSlider','Label','Max. BL',...
                  'SliderMinMaxVal',[0,20.0,TX(1).peakBLMax],...
                  'SliderStep',[0.005,0.020],'ValueFormat','%1.3f'};
UI(6).Callback = text2cell('%MaxBLCallback');

%% Additional UI for singleshot
% Fire singleshot
UI(7).Control = {'UserC2', 'Style', 'VsPushButton',...
                'Label', 'Acquire CF'};
UI(7).Callback = text2cell('%SingleshotActivateCallback');

% Add a control for wall filter thresholds
UI(8).Control = {'UserC1','Style','VsSlider','Label','SVD Thresh',...
                 'SliderMinMaxVal',[2, opts.singleshot_PRI-1, 15],...
                 'SliderStep',[0.01, 0.01],'ValueFormat','%.0f'};
UI(8).Callback = text2cell('%SVDSliderCallback');

% Let user specify filename
savefile = '';                              % Default to nothing in box
UI(9).Control = {'Style', 'edit',...        % Control to edit input text
                 'String', savefile,...     % Initial value
                 'Position', [165,220,275,22],...
                 'Tag','centerFreq',...
                 'BackgroundColor',[0.9,0.9,0.9],...
                 'Callback',{@filenameCallback}};
UI(9).Callback = {'filenameCallback.m',...
                  'function filenameCallback(hObject,eventdata)',...
                  ' ',...
                  'temp_filename = get(hObject,''String'');',...
                  'assignin(''base'',''savefile'',temp_filename);',...
                  'return'};
              
% Label filename box for clarity
UI(10).Control = {'Style', 'text',...
                  'String', 'Filename',...
                  'Position', [165,243,140,20]...
                  'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 9,...
                  'BackgroundColor', [0.8, 0.8, 0.8]};      

%% Specify external functions
EF(1).Function = text2cell('%SaveRFData');
EF(2).Function = text2cell('%SaveIQData');
EF(3).Function = text2cell('%ProcessDisplayDoppler');
EF(4).Function = text2cell('%BmodeFrameSave');
EF(5).Function = text2cell('%ResetStartingEvent');

%% Automatic VSX Execution:
% Uncomment the following line to automatically run VSX every time you run
% this SetUp script (note that if VSX finds the variable 'filename' in the
% Matlab workspace, it will load and run that file instead of prompting the
% user for the file to be used):

% filename = 'MUSIC_Colorflow';  VSX;

% Specify factor for converting sequenceRate to frameRate.
frameRateFactor = 3;

% Save all the structures to a .mat file.
save('MatFiles/MUSIC_Colorflow');
filename = 'MatFiles/MUSIC_Colorflow';

return


%% Define GUI callback and External functions
% **** Callback routines to be converted by text2cell function. ****

%% UI(1) - sensitivity cutoff
%SensCutOffCallback - Sensitivity cutoff change
ReconL = evalin('base', 'Recon');
for i = 1:size(ReconL,2)
    ReconL(i).senscutoff = UIValue;
end
assignin('base','Recon',ReconL);
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'Recon'};%SaveBmodeFrame
assignin('base','Control', Control);
return
%SensCutOffCallback

%% UI(2) - range change
%RangeChangeCallback - Range change
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.endDepthMm = UIValue;
        P.endDepth = UIValue*scaleToWvl;
    else
        P.endDepth = UIValue;
        P.endDepthMm = P.endDepth/scaleToWvl;
    end
end
assignin('base','P',P);

% No range change if in simulate mode 2.
if Resource.Parameters.simulateMode == 2
    return
end

% Modify PData and Regions for new depth.
PData = evalin('base','PData');
PData(1).Size(1) = ceil((P.endDepth-P.startDepth)/PData(1).PDelta(3)); % P.startDepth, P.endDepth and pdelta set PData.Size(1).
PData(1).Region = repmat(struct('Shape',struct( ...
                    'Name','Parallelogram',...
                    'Position',[0,0,P.startDepth],...
                    'width',5*Trans.spacing*Trans.numelements/(P.numRays),...
                    'height',P.endDepth-P.startDepth,...
                    'angle',0.0)),1,3*P.numRays);
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(P.numRays-1)):(63.5*Trans.spacing);
for n = 1:P.numRays, PData(1).Region(n).Shape.Position(1) = TxOrgX(n); end
m = P.numRays;
for n = 1:P.numRays
    if n<=8
        angle = -((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = -((P.numRays-n)/8)*P.dtheta;
    else
        angle = -P.dtheta;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
m = m + P.numRays;
% Define numRays steered right parallelogram regions, centered on TX beam origins. Adjust the angle
%   so that the steering goes to zero over 8 beams at the left and right edge.
for n = 1:P.numRays
    if n<=8
        angle = ((n-1)/8)*P.dtheta;
    elseif n>(P.numRays-8)
        angle = ((P.numRays-n)/8)*P.dtheta;
    else
        angle = P.dtheta;
    end
    PData(1).Region(n+m).Shape.Position(1) = TxOrgX(n);
    PData(1).Region(n+m).Shape.height = (P.endDepth-P.startDepth)/cos(angle);
    PData(1).Region(n+m).Shape.angle = angle;
end
PData(1).Region = computeRegions(PData(1));
assignin('base','PData',PData);
evalin('base','Resource.DisplayWindow(1).Position(4) = ceil(PData(1).Size(1)*PData(1).PDelta(3)/Resource.DisplayWindow(1).pdelta);');

% Update TXPD data of TX structures.
TX = evalin('base','TX');
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for ind = 1:steps
    TX(ind).TXPD = computeTXPD(TX(ind),PData(1));
    waitbar(ind/steps)
end
close(h)
assignin('base','TX',TX);

% Update Receive and TGC
Receive = evalin('base', 'Receive');
maxAcqLength = ceil(sqrt(P.endDepth^2 + ((Trans.numelements-1)*Trans.spacing)^2));
for i = 1:size(Receive,2)
    Receive(i).endDepth = maxAcqLength;
end
assignin('base','Receive',Receive);
evalin('base','TGC.rangeMax = P.endDepth;');
evalin('base','TGC.Waveform = computeTGCWaveform(TGC);');

Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'PData','InterBuffer','ImageBuffer','DisplayWindow','TX','Receive','TGC','Recon'};
assignin('base','Control', Control);
assignin('base', 'action', 'displayChange');
return
%RangeChangeCallback

%% UI(3) - transmit focus
%TxFocusCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.txFocus'));
    return
end
Trans = evalin('base','Trans');
Resource = evalin('base','Resource');
scaleToWvl = Trans.frequency/(Resource.Parameters.speedOfSound/1000);

P = evalin('base','P');
P.txFocus = UIValue;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm')
        P.txFocus = UIValue*scaleToWvl;
    end
end
assignin('base','P',P);%SaveBmodeFrame

% - Redefine event specific TX attributes for the new focus.
TX = evalin('base', 'TX');
PData = evalin('base','PData');
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*P.numRays;
for ind = 1:steps
    % write new focus value to TX
    TX(ind).focus = P.txFocus;
    TX(ind).Delay = computeTXDelays(TX(ind));
    TX(ind).TXPD = computeTXPD(TX(ind),PData(1));
    waitbar(ind/steps)
end
close(h)
assignin('base','TX',TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%TxFocusCallback

%% UI(4) - transmit aperture
%ApertureCallback
simMode = evalin('base','Resource.Parameters.simulateMode');
% No focus change if in simulate mode 2.
if simMode == 2
    set(hObject,'Value',evalin('base','P.numTx'));
    return
end
numTx = UIValue;
evalin('base',['P.numTx = ',int2str(numTx),';']);
Trans = evalin('base', 'Trans');
numRays = evalin('base', 'P.numRays');
dtheta = evalin('base','P.dtheta');
TX = evalin('base', 'TX');
PData = evalin('base','PData');
scaleToWvl = evalin('base','scaleToWvl');
TxOrgX = (-63.5*Trans.spacing):(127*Trans.spacing/(numRays-1)):(63.5*Trans.spacing);
% - Redefine event specific TX attributes for the new aperture.
h = waitbar(0,'Program TX parameters, please wait!');
steps = 3*numRays;
for n = 1:numRays  % specify P.numRays transmit events
    TX(n).Origin(1) = TxOrgX(n);
    % Compute transmit aperture apodization
    TX(n).Apod = +(((scaleToWvl*Trans.ElementPos(:,1))>(TxOrgX(n)-Trans.spacing*numTx/2))& ...
                 ((scaleToWvl*Trans.ElementPos(:,1))<(TxOrgX(n)+Trans.spacing*numTx/2)))';
    [RIndices,CIndices,V] = find(TX(n).Apod);
    V = kaiser(size(V,2),1);
    TX(n).Apod(CIndices) = V;
    % Compute transmit delays
    TX(n).Delay = computeTXDelays(TX(n));
    % Compute transmit pixel data
    TX(n).TXPD = computeTXPD(TX(n),PData(1));
    waitbar(n/steps);
end
m = numRays;
for n = 1:numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [-((n-1)/8)*dtheta,0.0];
    elseif n>(numRays-8)
        TX(n+m).Steer = [-((numRays-n)/8)*dtheta,0.0];
    else
        TX(n+m).Steer = [-dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Compute transmit pixel data
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData(1));
    waitbar((n+numRays)/steps);
end
m = m + numRays;
for n = 1:numRays
    TX(n+m).Origin(1) = TX(n).Origin(1);
    TX(n+m).Apod = TX(n).Apod;
    if n<=8
        TX(n+m).Steer = [((n-1)/8)*dtheta,0.0];
    elseif n>(numRays-8)
        TX(n+m).Steer = [((numRays-n)/8)*dtheta,0.0];
    else
        TX(n+m).Steer = [dtheta,0.0];
    end
    % Compute transmit delays
    TX(n+m).Delay = computeTXDelays(TX(n+m));
    % Compute transmit pixel data%SaveBmodeFrame
    TX(n+m).TXPD = computeTXPD(TX(n+m),PData(1));
    waitbar((n+2*numRays)/steps);
end
close(h)
assignin('base','TX', TX);
% Set Control command to update TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
return
%ApertureCallback

%% UI(5) - peak cutoff
%PeakCutOffCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakCutOff = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%PeakCutOffCallback

%% UI(6) - max burst length
%MaxBLCallback
TX = evalin('base', 'TX');
for i=1:size(TX,2)
    TX(i).peakBLMax = UIValue;
end
assignin('base','TX',TX);
% Set Control.Command to set TX
Control = evalin('base','Control');
Control.Command = 'update&Run';
Control.Parameters = {'TX'};
assignin('base','Control', Control);
%MaxBLCallback

%% UI(7) - Fire colorflow singleshot
%SingleshotActivateCallback
singleshot_start = evalin('base','singleshot_start');
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',singleshot_start};
assignin('base','Control', Control);
%SingleshotActivateCallback

%% UI(8) - Adjust SVD threshold
%SVDSliderCallback
opts = evalin('base', 'opts');
new_thresh = round(UIValue);
opts.filter_thresholds = [new_thresh opts.filter_thresholds(2)];
assignin('base','opts', opts);
%SVDSliderCallback

%% EF(1) - save RF if desired
%SaveRFData
save_raw_RF(RFData)
opts = evalin('base', 'opts');
savefile = evalin('base', 'savefile');

% Do nothing if not saving RF data
if ~opts.save_RF
    return
end

% Count acqs with matching filename, auto-increment
if strcmp(savefile, '')
    savefile = 'filename_empty';
end
flist = dir([opts.save_path '*' savefile '*RF.mat']);
fileCount = length(flist);

% If toggled on, grab RF data and save to specified location
savename = sprintf([savefile '_CF%02d_RF.mat'], fileCount);
filename = fullfile(opts.save_path, savename);
ticS = tic;
try
    savefast(filename, 'RFData');
catch
    fprintf('Please place "savefast" function from file exchange on Matlab path; using standard for now.\n')
    save(filename, 'RFData');
end
fprintf('\n%s saved: %.2f seconds.\n', savename, toc(ticS))
%SaveRFData

%% EF(2) - save IQ if desired
%SaveIQData
save_IQ_ensemble(IData, QData)
opts = evalin('base', 'opts');
savefile = evalin('base', 'savefile');

% Do nothing if not saving IQ data
if ~opts.save_IQ
    return
end

% Count acqs with matching filename, auto-increment
if strcmp(savefile, '')
    savefile = 'filename_empty';
end
flist = dir([opts.save_path '*' savefile '*IQ.mat']);
fileCount = length(flist);

% If toggled on, grab IQ data and save to specified location
savename = sprintf([savefile '_CF%02d_IQ.mat'], fileCount);
filename = fullfile(opts.save_path, savename);
ticS = tic;
IData = squeeze(IData); QData = squeeze(QData);
try
    savefast(filename, 'IData', 'QData');
catch
    fprintf('Please place "savefast" function from file exchange on Matlab path; using standard for now.\n')
    save(filename, 'IData', 'QData');
end
fprintf('%s saved: %.2f seconds.\n', savename, toc(ticS))
%SaveIQData

%% EF(3) - process and display Doppler results
%ProcessDisplayDoppler
process_Doppler(IData, QData)
ticS = tic;

% Get options
opts = evalin('base', 'opts');

% Identify where to save output; match IQ and RF core string if also saving
% raw data
savefile = evalin('base', 'savefile');
if strcmp(savefile, '')
    savefile = 'filename_empty';
end
flist = dir([opts.save_path '*' savefile '*IQ.mat']);
fileCount = length(flist);

% Setting up for saving pngs, figs; colorflow
dopPwrFigSave = sprintf([savefile '_CF%02d_Power.fig'], fileCount);
dopPwrPngSave = sprintf([savefile '_CF%02d_Power.png'], fileCount);
dopVelFigSave = sprintf([savefile '_CF%02d_Velocity.fig'], fileCount);
dopVelPngSave = sprintf([savefile '_CF%02d_Velocity.png'], fileCount);
dopRawSave = sprintf([savefile '_CF%02d_Doppler_output.mat'], fileCount);

% Prepare IQ Doppler processing
IQ = squeeze(complex(IData, QData));

% Wall filter
if ~strcmp(opts.filter_type, 'SVD')
    fprintf('Other wall filters not implemented yet, defaulting to SVD.\n')
else
    % Run SVD filter
    % Generate Casorati matrix from stack of IQ data
    sz = size(IQ);
    casorati = reshape(IQ, [sz(1)*sz(2),sz(3)]);
    
    % SVD filter to isolate signal of interest
    [V, D] = eig(cov(casorati));                         % Get eigenvals, vectors
    [~, ind] = sort(diag(D), 'descend'); vi = V(:,ind);  % Sort by descending eigenvalue
    lamu = casorati * vi;
    
    % Isolate desired signal
    startVector = opts.filter_thresholds(1); 
    endVector = opts.filter_thresholds(2);
    vi(:, 1:startVector-1) = 0; lamu(:, 1:startVector-1) = 0; % Reject first SVs, clutter
    if endVector < sz(3)
        vi(:, endVector:end) = 0; lamu(:, endVector:end) = 0; % Reject last SVs, noise
    end
    bloodSignal = lamu * ctranspose(vi);
    
    % Reshape back into image dimensions
    filtIQ = reshape(bloodSignal, [sz(1),sz(2),sz(3)]);
    sz = size(filtIQ);
end

% Get flow power
if strcmp(opts.power_processing, 'Cleaner')
    realSq = real(filtIQ).^2; 
    imagSq = imag(filtIQ).^2;
    sumIQ = realSq + imagSq; 
    pwr = sum(sumIQ, 3);
else
    pwr = abs(mean(filtIQ(:,:,1:end-1).*conj(filtIQ(:,:,2:end)),3));
end

% Get velocity
vel = angle(mean(filtIQ(:,:,1:end-1).*conj(filtIQ(:,:,2:end)),3)); % rad
vmax = opts.max_velocity;      % cm/sec
velCorr = -(vmax .* vel) / pi; % cm/sec

% Display flow power
pwrFig = figure('Units', 'normalized', 'Position', [0.37, 0.25, 0.25, 0.5], ...
                'Color', 'k', 'Name', dopPwrFigSave);
imagesc(10*log10(pwr)); colormap gray; axis image; 
h = colorbar('Color', 'w'); ylabel(h, 'dB', 'Color', 'w'); caxis(opts.power_display);
title('Flow Power', 'Color', 'w');
set(gca,'XColor','w','YColor','w');
xtx = get(gca, 'xtick'); ztx = get(gca, 'ytick');
mmx = opts.mmPerPxX_dop .* xtx; mmz = opts.mmPerPxZ_dop .* ztx;
set(gca,'xticklabel',split(sprintf('%.1f|',mmx), '|'))
set(gca,'yticklabel',split(sprintf('%.1f|',mmz), '|'))
xlabel('mm', 'Color', 'w');ylabel('mm', 'Color', 'w')
set(gcf, 'InvertHardcopy', 'off');

% Define vessel mask to threshold out noise
vesselMsk = sigmf(10*log10(pwr), [20, opts.velocity_power_threshold]);  % For thresholding velocity img

% Display velocity
velFig = figure('Units', 'normalized', 'Position', [0.1, 0.25, 0.25, 0.5], ...
                'Color', 'k', 'Name', dopVelFigSave);
load('CFI-colormap.mat')
imagesc(medfilt2(velCorr.*imbinarize(vesselMsk))); colormap(CFI); axis image; 
h = colorbar('Color', 'w'); ylabel(h, 'cm/sec', 'Color', 'w'); caxis(opts.velocity_display);
title('Flow Velocity', 'Color', 'w');
set(gca,'XColor','w','YColor','w');
xtx = get(gca, 'xtick'); ztx = get(gca, 'ytick');
mmx = opts.mmPerPxX_dop .* xtx; mmz = opts.mmPerPxZ_dop .* ztx;
set(gca,'xticklabel',split(sprintf('%.1f|',mmx), '|'))
set(gca,'yticklabel',split(sprintf('%.1f|',mmz), '|'))
xlabel('mm', 'Color', 'w');ylabel('mm', 'Color', 'w')
set(gcf, 'InvertHardcopy', 'off');

% Save raw data if desired
if opts.save_Doppler_output
    try
        savefast(fullfile(opts.save_path, dopRawSave), 'pwr', 'vel', 'velCorr', 'opts');
    catch
        save(fullfile(opts.save_path, dopRawSave), 'pwr', 'vel', 'velCorr', 'opts');
    end
    fprintf('Doppler output saved to %s.\n', dopRawSave);
end

% Save images if desired
if opts.save_images
    saveas(pwrFig, fullfile(opts.save_path, dopPwrFigSave));
    saveas(pwrFig, fullfile(opts.save_path, dopPwrPngSave));
    saveas(velFig, fullfile(opts.save_path, dopVelFigSave));
    saveas(velFig, fullfile(opts.save_path, dopVelPngSave));
    fprintf('Doppler images saved.\n');
end

fprintf('Doppler processing complete: %.2f seconds.\n', toc(ticS))
%ProcessDisplayDoppler

%% EF(4) - save Bmode frame associated with CF acquisition
%BmodeFrameSave
save_Bmode(imgData)
ticS = tic;

% Get options
opts = evalin('base', 'opts');
PData = evalin('base', 'PData');

% Identify where to save output; match IQ and RF core string if also saving
% raw data
savefile = evalin('base', 'savefile');
if strcmp(savefile, '')
    savefile = 'filename_empty';
end
flist = dir([opts.save_path '*' savefile '*IQ.mat']);
fileCount = length(flist);

bmodeIntSave = sprintf([savefile '_CF%02d_Bmode_raw_intensity.mat'], fileCount);
bmodeFigSave = sprintf([savefile '_CF%02d_Bmode.fig'], fileCount);
bmodePngSave = sprintf([savefile '_CF%02d_Bmode.png'], fileCount);

% Display Bmode image for saving
imgData = squeeze(imgData);
lastFrame = imgData(:,:,end-3:end);
comprImg = 20*log10(lastFrame(:,:,end));
maxInt = max(comprImg(:));
scal = [maxInt-opts.bmode_display_dynamic_range, maxInt]; % 60 dB display dynamic range by default
figure('Units', 'normalized', 'Position', [0.64, 0.25, 0.25, 0.5], ...
       'Color', 'k', 'Name', bmodeFigSave); 
imagesc(comprImg); colormap gray; axis image; 
h = colorbar('Color', 'w'); ylabel(h, 'dB', 'Color', 'w'); caxis(scal);
title('Bmode, log compressed', 'Color', 'w');
set(gca,'XColor','w','YColor','w');
xtx = get(gca, 'xtick'); ztx = get(gca, 'ytick');
mmx = opts.mmPerPxX .* xtx; mmz = opts.mmPerPxZ .* ztx;
set(gca,'xticklabel',split(sprintf('%.1f|',mmx), '|'))
set(gca,'yticklabel',split(sprintf('%.1f|',mmz), '|'))
xlabel('mm', 'Color', 'w');ylabel('mm', 'Color', 'w')
set(gcf, 'InvertHardcopy', 'off');

% Exit if not saving
if ~opts.save_bmode_images
    return
end

% Otherwise, save
try
    savefast(fullfile(opts.save_path, bmodeIntSave), 'imgData');
catch
    save(fullfile(opts.save_path, bmodeIntSave), 'imgData');
end
saveas(gcf, fullfile(opts.save_path, bmodeFigSave));
saveas(gcf, fullfile(opts.save_path, bmodePngSave));

fprintf('%s saved: %.2f seconds.\n\n', bmodePngSave, toc(ticS))
%BmodeFrameSave

%% EF(5) - reset start event to jump back to Bmode after singleshot
%ResetStartingEvent
reset_start_event()
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Parameters',1,'startEvent',1};
assignin('base','Control', Control);
%ResetStartingEvent


