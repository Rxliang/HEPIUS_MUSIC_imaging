function MSet = readrfbwtfn(inFile, MSpec)

tSel_s = MSpec.tSel_s;
lslrt(inFile)
fid = fopen(inFile, 'r');

MSet.fileVersion = fread(fid, 1, 'uint8');
MSet.numLinesTotal = fread(fid, 1, 'uint16');
MSet.numEnsemble = fread(fid, 1, 'uint16');
MSet.numLines = MSet.numLinesTotal/MSet.numEnsemble
MSet.numLinesPerFrame = fread(fid, 1, 'uint16');
MSet.freq_MHz = fread(fid, 1, 'single');
MSet.frameRate_Hz = fread(fid, 1, 'single');
MSet.lumenLateral_mm = fread(fid, 1, 'single');
MSet.lumenDepth_mm = fread(fid, 1, 'single');
MSet.lumenWidth_mm = fread(fid, 1, 'single');
MSet.lumenHeight_mm = fread(fid, 1, 'single');
MSet.lenLine = fread(fid, 1, 'uint16');
MSet.zAx_mm = fread(fid, MSet.lenLine, 'single');
MSet.dateNumVec = fread(fid, MSet.numLines, 'double');
MRaw = fread(fid, Inf, 'int16');

fclose(fid);

MSet.xGate_mm = MSet.lumenLateral_mm;
MSet.zGate_mm = MSet.lumenDepth_mm+MSet.lumenHeight_mm/2;
lenRaw = length(MRaw);
%numLines = floor(lenRaw/MSet.lenLine);

MPre = reshape(MRaw, MSet.lenLine, MSet.numLinesTotal);

[atten, comp] = usattenfn(MSet.zAx_mm, MSet.freq_MHz, MSpec.tissueAtten_dBpMHz);
COMP = repmat(comp(:), 1, size(MPre,2));

MSet.M = MPre .* COMP;

MSet.lineRate_Hz = MSet.frameRate_Hz*MSet.numEnsemble;
MSet.tAx_s = 0:1/MSet.lineRate_Hz:(MSet.numLinesTotal-1)/MSet.lineRate_Hz;
MSet.MSpec = MSpec;
MSet.inFile = inFile;
MSet.figHandle = [];
MSet.imgHandle = [];
MSet.lineHandle = [];


dnVecAug = [MSet.dateNumVec;  MSet.dateNumVec(end)+ ...
            mean(diff(MSet.dateNumVec))];

dnVecInterp = [];
    
for r = 1:length(dnVecAug)-1
  intVec = linspace(dnVecAug(r), dnVecAug(r+1), MSet.numEnsemble+1);
  dnVecInterp = [dnVecInterp intVec(1:end-1)];
end

MSet.dateNumVecInterp = dnVecInterp;

if isfield(MSpec, 'PAug')
  MSet.pressureInterp1_mmHg = interp1(MSpec.PAug.datenum, MSpec.PAug.p1, ...
                                         MSet.dateNumVecInterp);
  MSet.pressureInterp2_mmHg = interp1(MSpec.PAug.datenum, MSpec.PAug.p2, ...
                                      MSet.dateNumVecInterp);
    
  MSet.meanPressure_mmHg = mean([MSet.pressureInterp1_mmHg;
                                 MSet.pressureInterp2_mmHg],1);
                      
end

if MSpec.plt
  MSet.figHandle = figure(MSpec.plt)
  clf
  colormap(parula)
  MSet.imgHandle = imagesc(MSet.tAx_s, MSet.zAx_mm-MSpec.startDepthCorr_mm, ...
                           log(sqrt(abs(MSet.M))))
  yl = ylim;
  xlabel('time (s)');
  ylabel('depth (mm)');
  hold on
  MSet.lineHandle = line(xlim, MSet.zGate_mm*[1 1]-MSpec.startDepthCorr_mm, 'linestyle',['-' ...
                      '-'],'color', 'w', 'linewidth',2);

  if isfield(MSpec, 'PAug')
    yyaxis right;
    plot(MSet.tAx_s, ...
         MSet.meanPressure_mmHg, 'w', 'linewidth', 2);
    ylabel('pressure (mmHg)');
  end
  
  hold off
  
end
