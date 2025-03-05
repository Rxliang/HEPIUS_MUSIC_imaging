function [Process, proc, Resource, evParms] = sdcprocfn(evParms, Resource)

PDataIndB = evParms.ind.PDataIndB;
doCFI = evParms.flag.doCFI;
numSD = evParms.gate.numSD;
interBufferIndSDStart = evParms.ind.interBufferIndSDStart;
interBufferIndCFI = evParms.ind.interBufferIndCFI;
imageBufferIndCFI = evParms.ind.imageBufferIndCFI;
imageBufferIndSDLarge = evParms.ind.imageBufferIndSDLarge;
imageBufferIndSDLargeIntermed = evParms.ind.imageBufferIndSDLargeIntermed;
PIndCFI = evParms.ind.PIndCFI;
interBufferIndSDLarge = evParms.ind.interBufferIndSDLarge;
interBufferIndSDLargeOut = evParms.ind.interBufferIndSDLargeOut;

% Specify Process structure arrays.
proc = [];
proc.procIndB = 1;
persistenceB = 20;
Process = [];
Process(proc.procIndB).classname = 'Image';
Process(proc.procIndB).method = 'imageDisplay';
% if we are not waiting for CFI overlay, show image now. Note, if
% displayBNow is a logical type, does not enable display. learned
% the very hard way

displayBNow= double(~doCFI);

if evParms.state.SDOnlyInPlace
  % reduce image gain for SDOnlyInPlace because excitation have
  % more half-cycles
  evParms.ev.BImGain = evParms.largeSDParms.BImGainScale;
else
  evParms.ev.BImGain = 1;
end

if evParms.flag.BLineWT
  evParms.ev.BImGain = 0.1 + ...
      0.75*(evParms.TW(evParms.ind.TWIndB).Parameters(1)/evParms.Trans.frequency-1);

  systemID = getenv('BNS_SYSTEM_ID');
  if strcmp(systemID, 'vera64')
    evParms.ev.BImGain = evParms.ev.BImGain*3;
  end
  
end

% this to be used along with CFI
Process(proc.procIndB).Parameters = {'imgbufnum', 1, ... % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum', PDataIndB,...
                         'pgain', evParms.ev.BImGain, ...
                         'reject', 2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel', persistenceB,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod','lowerHalf',... 
                         'display', displayBNow, ... % display image after processing
                         'displayWindow',1};
% NOTE: mappingMethod lowerHalf used with B/CFI, not B-only
proc.procIndSDExtStart = 2;

for i = 1:numSD
  Process(proc.procIndSDExtStart+i-1).classname = 'External';
  timingOnly=0;
  if ~timingOnly
    Process(proc.procIndSDExtStart+i-1).method = [evParms.SD.functionName '_' ...
                        num2str(i)];
    Process(proc.procIndSDExtStart+i-1).Parameters = {'srcbuffer','inter', ... 
                        'srcbufnum', interBufferIndSDStart+i-1, ...
                        'srcframenum',1, ...
                        'dstbuffer','none'};
  else
    uscmsg(evParms.state.execState);
    % development code
    % timing function for measuring SD invocation intervals
    Process(proc.procIndSDExtStart+i-1).method = ['extclockfn_' num2str(i)];
    Process(proc.procIndSDExtStart+i-1).Parameters = [];
  end
end

if evParms.flag.largeSD 
  proc.procIndSDLarge = proc.procIndSDExtStart+numSD;
  numSDProc = numSD+1;  
  if evParms.largeSDParms.doSDLargeInternalProc
    % process structure for 1st CFI Doppler ensemble
    Process(proc.procIndSDLarge).classname = 'Doppler';
    switch evParms.largeSDParms.SDLargeCDIMode
      case 'freq'
        Process(proc.procIndSDLarge).method = 'computeCFIFreqEst';
      case 'power'
        Process(proc.procIndSDLarge).method = 'computeCFIPowerEst';
      otherwise
        error('Unsupported SD overlay CDI state');
    end
      
    if evParms.largeSDParms.procSDLargeOut
      imgBufDest = [imageBufferIndSDLargeIntermed,1];
    else
      imgBufDest = [imageBufferIndSDLarge,-2];
    end
      
    if evParms.state.SDOnlyInPlace
      % SDLarge keeps using the 2-spaced PRIs, so it has half the
      % SD freq
      procSDLargePRF = evParms.ev.dopPRF/2;
    else
      procSDLargePRF = evParms.ev.dopPRF/evParms.largeSDParms.priSkip;
    end
      
    numPages = evParms.largeSDParms.endPRI-evParms.largeSDParms.startPRI+1;
    Process(proc.procIndSDLarge).Parameters = ...
        {'IntBufSrc',[interBufferIndSDLarge,1],...
         'SrcPages',[evParms.largeSDParms.startPRI, numPages],... 
         'ImgBufDest', imgBufDest,...                     
         'pdatanum', evParms.gate.PDataIndSDLarge,...    
         'prf', procSDLargePRF,...  % Doppler PRF in Hz
         'pwrThreshold', evParms.largeSDParms.pwrThres,...
         'wallFilter', evParms.largeSDParms.wallFilter, ...
         'maxPower', evParms.largeSDParms.maxPower,...
         'postFilter',0};
    %                            'IntBufDest',
    %                            [interBufferIndSDLargeOut,1],...
    %                            'ImgBufDest', [imageBufferIndSDLarge,-2], ... 

  else
    % development code for external processing
    if 0
       uscmsg(execState);
       % external processing  
       Process(proc.procIndSDLarge).classname = 'External';
       Process(proc.procIndSDLarge).method = 'cdisdlprocessing_jon';
       Process(proc.procIndSDLarge).Parameters = {'srcbuffer','inter',...  % name of buffer to process.
                           'srcbufnum',interBufferIndSDLarge,...
                           'srcframenum',1,...
                           'dstbuffer','image' ...
                           'dstbufnum',imageBufferIndSDLargeIntermed, ...
                           'dstframenum',1};
    end  
    
    if 0
       uscmsg(execState);
       Process(proc.procIndSDLarge).classname = 'External';
       Process(proc.procIndSDLarge).method = ['sdlargefn'];
       Process(proc.procIndSDLarge).Parameters = {'srcbuffer','inter', ... % name of buffer to process.
                           'srcbufnum', interBufferIndSDStart+numSD, ...
                           'srcframenum',1, ...
                           'dstbuffer','image', ...
                           'dstbufnum',imageBufferIndSDLargeIntermed, ...
                           'dstframenum', 1}; % was -1why was it
                                              % -2? 
    end
  end
  
else % no SDLarge
  numSDProc = numSD;
end

% adapted from widebeamDoppler
proc.procIndCFI = 1+numSDProc+1;
  
if evParms.state.CDIState==1 | evParms.state.CDIState==3
    Resource.DisplayWindow(1).Colormap = grayscaleCFImap;
    imSrcData = 'signedColor';
    imPersistMethod = 'dynamic';
    if evParms.state.CDIState==1
      imPersistLevel = evParms.P(PIndCFI).persfreq;
    else
      imPersistLevel = evParms.largeSDParms.persfreq;
    end
    
end

if evParms.state.CDIState==2 | evParms.state.CDIState==4
  Resource.DisplayWindow(1).Colormap = grayscaleCPAmap;
  imSrcData = 'unsignedColor';
  imPersistMethod = 'simple';
  if evParms.state.CDIState==2
    imPersistLevel = evParms.P(PIndCFI).perspwr;
  else
    imPersistLevel = evParms.largeSDParms.perspwr;
  end
end

if evParms.flag.doCFIInternalProc
  % process structure for 1st CFI Doppler ensemble
  Process(proc.procIndCFI).classname = 'Doppler';
  
  if evParms.state.CDIState==1
    Process(proc.procIndCFI).method = 'computeCFIFreqEst';
  else
    Process(proc.procIndCFI).method = 'computeCFIPowerEst';
  end
  
  % number of (Inter?)buffer to process:
  % start frame number in source buffer:
  Process(proc.procIndCFI).Parameters = {'IntBufSrc',[interBufferIndCFI,1],...
                      'SrcPages',[3,evParms.P(PIndCFI).dopPRIs-2],...   
                      'ImgBufDest', [imageBufferIndCFI,-1],...
                      'pdatanum', evParms.ind.PDataIndCFI1,...    
                      'prf', evParms.P(PIndCFI).dopPRF,...  % Doppler PRF in Hz
                      'wallFilter', evParms.CFI.wallFilter,...  % 1 -> quadratic regression
                      'pwrThreshold',evParms.P(PIndCFI).pwrThres,...
                      'maxPower',50,...
                      'postFilter',1};
else
  uscmsg(evParms.state.execState);
  Process(proc.procIndCFI).classname = 'External';
  Process(proc.procIndCFI).method = 'cdiprocessing_jon';
  Process(proc.procIndCFI).Parameters = {'srcbuffer','inter',...  
                      'srcbufnum',interBufferIndCFI,...
                      'srcframenum',1,...
                      'dstbuffer','image' ...
                      'dstbufnum',imageBufferIndCFI, ...
                      'dstframenum',-2};
  % Both options write to the last frame in the destination buffer, but -2 will increment the last frame pointer and
  % the Receive.ImageBuffer.lastFrame attribute before writing
end

proc.procIndCFIIm = proc.procIndCFI+1;
Process(proc.procIndCFIIm).classname = 'Image';
Process(proc.procIndCFIIm).method = 'imageDisplay';
Process(proc.procIndCFIIm).Parameters = ...
    {'imgbufnum',imageBufferIndCFI,...   % number of buffer to process.
     'framenum',-1,...   % (-1 => lastFrame)
     'srcData', imSrcData,... % type of data to display.
     'pdatanum',evParms.ind.PDataIndCFI1,...    % number of PData structure to use
     'pgain',1,...            % pgain is image processing gain
     'reject',2,...      % reject level
     'persistMethod',imPersistMethod, ...
     'persistLevel', imPersistLevel, ...
     'interpMethod','4pt',...  %method of interp. (1=4pt)
     'grainRemoval','none',...
     'processMethod', evParms.CFI.processMethod,...
     'averageMethod','none',...
     'compressMethod','power',...
     'compressFactor', evParms.CFI.compressFactor,...
     'mappingMethod','upperHalf',...
     'threshold', evParms.P(PIndCFI).cpl, ...
     'display',1,...      % display image after processing
     'displayWindow',1};

% EF1 is external function for ROI plot
proc.procIndCFIROIPlot = proc.procIndCFIIm+1;
Process(proc.procIndCFIROIPlot).classname = 'External';
Process(proc.procIndCFIROIPlot).method = 'nicpseqwt_roiplot';
Process(proc.procIndCFIROIPlot).Parameters = {'srcbuffer','none'};

if 1 % always make process % doCFI
  proc.procIndBNow = 6+numSDProc-1;
else
  proc.procIndBNow = 1+numSDProc+1;
end
  
Process(proc.procIndBNow).classname = 'Image';
Process(proc.procIndBNow).method = 'imageDisplay';
Process(proc.procIndBNow).Parameters = {'imgbufnum', 1, ... % number of buffer to process.
                         'framenum',-1,...   % (-1 => lastFrame)
                         'pdatanum', PDataIndB,...    % number of PData structure to use
                         'pgain', evParms.ev.BImGain,... % pgain is image processing gain
                         'reject', 2,...      % reject level 
                         'persistMethod','simple',...
                         'persistLevel', persistenceB,...
                         'interpMethod','4pt',...  %method of interp. (1=4pt)
                         'grainRemoval','none',...
                         'processMethod','none',...
                         'averageMethod','none',...
                         'compressMethod','power',...
                         'compressFactor',40,...
                         'mappingMethod', 'lowerHalf', ... % trying
                         'display', 1, ... % display image after processing
                         'displayWindow',1};

proc.procIndSDLargeIm = proc.procIndBNow+1;

Process(proc.procIndSDLargeIm).classname = 'Image';
Process(proc.procIndSDLargeIm).method = 'imageDisplay';
Process(proc.procIndSDLargeIm).Parameters = ...
    {'imgbufnum',imageBufferIndSDLarge,...   % number of buffer to process.
     'framenum',-1,...   % (-1 => lastFrame)
     'srcData', imSrcData,... % type of data to display.
     'pdatanum', evParms.gate.PDataIndSDLarge,...    % number of PData structure to use
     'pgain',1,...            % pgain is image processing gain
     'reject',50,...      % reject level, percentage of lower, was 2
     'persistMethod', imPersistMethod, ...
     'persistLevel', imPersistLevel, ...
     'interpMethod','4pt',...  %method of interp. (1=4pt)
     'grainRemoval','none',...
     'processMethod','none',...
     'averageMethod','none', ... %'runAverage3',...
     'compressMethod','power',...
     'compressFactor', 20,... % 40 is sqrt, 20 is linear
     'mappingMethod','upperHalf',... %     'mappingMethod','upperHalf',...
     'threshold', evParms.largeSDParms.cpl,...
     'display',1,...      % display image after processing
     'displayWindow',1};

procCnt = proc.procIndSDLargeIm;

procCnt = procCnt+1;
proc.procIndSDLargeROIPlot = procCnt;;
Process(proc.procIndSDLargeROIPlot).classname = 'External';
Process(proc.procIndSDLargeROIPlot).method = 'ROIplot';
Process(proc.procIndSDLargeROIPlot).Parameters = {'srcbuffer','none'};

procCnt = procCnt+1;
proc.procIndSDLargeOut = procCnt;

Process(proc.procIndSDLargeOut).classname = 'External';
  Process(proc.procIndSDLargeOut).method = ['sdlargeoutfn'];
  Process(proc.procIndSDLargeOut).Parameters = ...
      {'srcbuffer','image', ... 
       'srcbufnum', imageBufferIndSDLargeIntermed, ...
       'srcframenum', 1, ...
       'dstbuffer', 'image',...
       'dstbufnum', imageBufferIndSDLarge, ...
       'dstframenum',-2};  


procCnt = procCnt+1;
proc.procIndSDLargeZeroOut = procCnt;

% return zeros for SDLarge
Process(proc.procIndSDLargeZeroOut).classname = 'External';
Process(proc.procIndSDLargeZeroOut).method = ['sdlargeoutzerofn'];
Process(proc.procIndSDLargeZeroOut).Parameters = ...
    {'srcbuffer','image', ... 
     'srcbufnum', imageBufferIndSDLargeIntermed, ...
     'srcframenum', 1, ...
     'dstbuffer', 'image',...
     'dstbufnum', imageBufferIndSDLarge, ...
     'dstframenum',-2 };  

% development code
if evParms.largeSDParms.largeSDSave;
    %  uscmsg(execState);
  procCnt = procCnt+1;
  proc.procIndSDLargeSave = procCnt;

  Process(proc.procIndSDLargeSave).classname = 'External';
  Process(proc.procIndSDLargeSave).method = ['sdlargesaveiqfn'];
  Process(proc.procIndSDLargeSave).Parameters = ...
      {'srcbuffer','inter', ... 
       'srcbufnum', interBufferIndSDLarge, ...
       'srcframenum', 1, ...
       'dstbuffer', 'none'};
end

      
% development code
if 0
  uscmsg(execState);
  Process(proc.procIndSDLargeOut).classname = 'External';
  Process(proc.procIndSDLargeOut).method = ['sdlargeoutfn'];
  Process(proc.procIndSDLargeOut).Parameters = ...
      {'srcbuffer','inter', ... 
       'srcbufnum', interBufferIndSDLargeOut, ...
       'srcframenum',1, ...
       'dstbuffer','none'};
end

% development code
if evParms.flag.BLineWT
  uscmsg(evParms.state.execState);
  procCnt = procCnt+1;
  proc.procIndBWT = procCnt;
  
  Process(proc.procIndBWT).classname = 'External';
  Process(proc.procIndBWT).method = 'procrflinefn';
  Process(proc.procIndBWT).Parameters = ...
                          {'srcbuffer',...
                           'receive',... % name of buffer to process.
                           'srcbufnum',1,... 
                           'srcframenum',-1,... 
                           'dstbuffer','none'};
end

% development code
if evParms.flag.wtCapability
  uscmsg(evParms.state.execState);
  procCnt = procCnt+1;
  proc.procIndWT = procCnt;
  
  Process(proc.procIndWT).classname = 'External';
  Process(proc.procIndWT).method = 'procwtfn';
  Process(proc.procIndWT).Parameters = ...
                          {'srcbuffer',...
                           'inter',... 
                           'srcbufnum', evParms.ind.interBufferIndWT,... 
                           'srcframenum', 1,... 
                           'dstbuffer','none'};
     
end

if evParms.flag.canSaveImage
  procCnt = procCnt+1;
  proc.procIndImSave = procCnt;
  
  Process(proc.procIndImSave).classname = 'External';
  %  Process(proc.procIndImSave).method = 'saveimwrapfn';
  Process(proc.procIndImSave).method = 'saveimfn';
  Process(proc.procIndImSave).Parameters = ...
                          {'srcbuffer',...
                           'imageP',... 
                           'srcbufnum', evParms.ind.imageBufferIndB, ...
                           'srcframenum', -1,... 
                           'dstbuffer','none'};

  %                           'srcbufnum', evParms.ind.imageBufferIndSDLarge, ...
end


end

