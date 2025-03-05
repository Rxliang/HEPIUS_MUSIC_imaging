% Version info
version = versionfn;
 
% default settings
% initial voltage [B/SD CFI] 
VInit = [1.6 1.6];
%VInit = [0 0];
% maximum voltage allowed [B/SD CFI] 
VMax = [96 96];


hn = gethostnamefn;

% to test on 256 channel machine
systemID = getenv('BNS_SYSTEM_ID');
if strcmp(systemID, 'vera256')
    disp('*** Settings for 256 channel machine applied ***');
    UTA260MUX=0;
    usem5_64=0;
end

if strcmp(systemID, 'veraBCM')
    disp('*** Settings for 64 channel machine, UTA-260S applied ***');
    % was 20x15
    SDLarge.maxSpgArea_mm2 = 30*20; % conservative, overwrites default
    UTA260MUX=1;
    usem5_64=1;
end

if  strcmp(hn, 'FULL-VVRILHIUJE')
    disp('*** Settings for 64 channel machine, UTA-260MUX=1 applied ***');
    % was 20x15
    SDLarge.maxSpgArea_mm2 = 30*20; % conservative, overwrites default
                                    %UTA260MUX=0;
    usem5_64=0;
    usem5=0;
    useMusic=1;
    UTA260MUX=2;
end

if  strcmp(hn,'FULL-HFMK7TQI8O')
    disp('*** Settings for 256 channel machine non-GE, UTA-260MUX=1 applied ***');
    % was 20x15
    SDLarge.maxSpgArea_mm2 = 30*20; % conservative, overwrites default
                                    %UTA260MUX=0;
    usem5_64=0;
    usem5=0;
    useMusic=3;
    UTA260MUX=3;
end

seqState = 'BSD'; % alternative is BCDI
SDOnlyInPlace=0;

disp('*** settings for sdlarge development with M5ScD');
%  PRIAcq_uSec = 111; % highest needed for B/SD duplex
%  PRIAcq_uSec = 150.5; % highest possible for approx 10cm depth 
PRIAcq_uSec = 221; % 4.5kHz with one wave type
PRIAcqInPlace_uSec = 221;

%PRIAcq_uSec = 300; % 4.5kHz with one wave type
%PRIAcqInPlace_uSec = 300;


forceRecord=0; % better to use semaphore to control this, setsemaon/off
largeSD=1;

largeSDParms.largeSDSave=0;
largeSDParms.doSDLargeInternalProc=1;
largeSDParms.SDLargeCDIMode =  'power';
largeSDParms.startPRI = 3;   
largeSDParms.endPRI = 79; % 80 or above will crash      
largeSDParms.wallFilter = 'FIRHigh';

if demoStableCase1 
    largeSDParms.startPRI = 3;
    largeSDParms.endPRI = 79; % 80 or above will crash      
    largeSDParms.wallFilter = 'FIRHigh';
end

if demoStableCase2
    largeSDParms.startPRI = 3;
    largeSDParms.endPRI = 100; % 80 or above will crash      
    largeSDParms.wallFilter = 'none';
end

if demoUnstableCase
    largeSDParms.startPRI = 3;
    largeSDParms.endPRI = 80; % 80 or above will crash      
    largeSDParms.wallFilter = 'FIRHigh';
end

largeSDParms.pwrThres = 0.035;
largeSDParms.procSDLargeOut=1; % process IQ output of computeCFIEst
largeSDParms.SDLargeDisplayOverlay=1;
% B-overwrite CDI pxi val thresh out of max 255
largeSDParms.cpl = 100; 
largeSDParms.maxPower=100; % JSM was 10
largeSDParms.perspwr=40; % persistence for power Doppler
largeSDParms.persfreq=40; % persistence for velocity Doppler
largeSDParms.priSkip=1; % enable us to skip pri's for overlay

dopWFCutoffNorm = 0.025*2; % SD wall filter cutoff as frac of PRF

SDLarge.zSDStart_mm = 0;
SDLarge.zSDHeight_mm = 85;
SDLarge.zSDEnd_mm = SDLarge.zSDStart_mm+SDLarge.zSDHeight_mm;

forceHighPRFModeForRecordingIQ = 1;

if 0
rcvParms.lnaGain_dB.B = 18;
rcvParms.lnaZinSel.B = 31;

rcvParms.lnaGain_dB.BWT = 24;
rcvParms.lnaZinSel.BWT = 31;

rcvParms.lnaGain_dB.CFI = 24;
rcvParms.lnaZinSel.CFI = 31;

rcvParms.lnaGain_dB.SD = 24;
rcvParms.lnaZinSel.SD = 31;
end


if 1
rcvParms.lnaGain_dB.B = 15;
rcvParms.lnaZinSel.B = 0;

rcvParms.lnaGain_dB.BWT = 15;
rcvParms.lnaZinSel.BWT = 0;

rcvParms.lnaGain_dB.CFI = 15;
rcvParms.lnaZinSel.CFI = 0;

rcvParms.lnaGain_dB.SD = 15;
rcvParms.lnaZinSel.SD = 0;
end

