mfile = mfilename;

vtData = getenv('BNS_DATA');

tuneParms.maxDia_mm = 3.5;
tuneParms.maxDiaGrowFactor = 1;
tuneParms.aMin_mm = tuneParms.maxDia_mm/4;
tuneParms.aMax_mm = tuneParms.maxDia_mm/2;
tuneParms.bMin_mm = tuneParms.maxDia_mm/4;
tuneParms.bMax_mm = tuneParms.maxDia_mm/2;
tuneParms.r0DeltaAsFracMaxDia = 4; 
tuneParms.imValMin = 0;
tuneParms.imValMax = 255;
plt = [];
plt.flag=0;
plt.figNo=2;
plt.solve=1;

srcScript = 'widebeamdoppler_l8';
dataPath = fullfile(vtData, srcScript);

h5FilePre = [srcScript '_'];
h5FilePost = '20200612_1555';
h5FilePost = '20200613_2351';
h5FilePost = '20200616_1728';
h5FilePost = '20200616_1753';
h5FilePost = '20200627_0103'; % train  frames = [112:445];
h5FilePost = '20200629_2345'; frames = [40:858]; % train
h5FilePost = '20200630_0004'; frames = [1:858]; % train


h5FileCDIPre = [srcScript '_'];
h5FileCDIPost = ['cdi_' h5FilePost];

h5File = fullfile(dataPath, [h5FilePre h5FilePost '.h5']);
lslrt(h5File);

h5FileCDI = fullfile(dataPath, [h5FileCDIPre h5FileCDIPost '.h5']);
lslrt(h5FileCDI);

[pth, fln, ext] = fileparts(h5File);

%outFile = fullfile(pth, [ mfile '_' fln '.mat']);

BFrame = h5read(h5File, '/data/B');
CDIFrame = h5read(h5FileCDI, '/data/CDI');
datenumFlow = h5read(h5FileCDI, '/data/datenumFrame');
tPre_dn = datenumFlow-datenumFlow(1);
t_s = tPre_dn*86400;

% interpolate to highest isotropic res
lambda_mm = h5read(h5File, '/header/lambda_mm');
d_mm = 0.5 * lambda_mm;

colMap = grayscaleCPAmap;

flowImSer = [];
fifoLen = 4;
cnt = 0;
cntFIFO = 1;
out = [];
out.flowImSto = [];
out.BImSto = [];
out.flowSum = [];
out.t_s = [];
out.frames = frames;
for q = 82:size(BFrame,3) %400:735 %
  q
  showImPre = BFrame(:,:,q).^0.6;
  BImStruct = readimh5fn(h5File, showImPre, d_mm);
  showIm = BImStruct.imInterp;
  showImScale = (showIm-minall(showIm))/(maxall(showIm)-minall(showIm))*127;
  flowImPre = CDIFrame(:,:,q);
  tThis_s = t_s(q);
  flowImStruct = readimh5fn(h5FileCDI, flowImPre, d_mm, ...
                            BImStruct.XI_mm, BImStruct.ZI_mm);
  
  flowImPre =  flowImStruct.imInterp;

  pos = mod(cnt, fifoLen)+1;
  flowImSer(:,:, pos) = flowImPre;
%  flowIm = mean(flowImSer(:,:,1:cntFIFO),3);
  flowIm = max(flowImSer(:,:,1:cntFIFO),[],3);
  cntFIFO = min(cntFIFO+1, fifoLen);
  cnt = cnt+1;
  ind = find(flowIm < showImScale);
  flowIm(ind)=0;
  ind = find(flowIm);
  flowIm(ind) = flowIm(ind)+127;

  out.flowImSto(:,:,cnt) = flowIm;
  out.BImSto(:,:,cnt) = showImScale;
  out.flowSum(cnt)= sumall(flowIm);
  out.t_s(cnt)= tThis_s;
  duplexIm = showImScale + flowIm;
  
  figure(1)
  imshow(duplexIm, [0 256], 'colorMap', colMap, 'xdata', ...
         BImStruct.xAxInterp_mm, 'ydata', BImStruct.zAxInterp_mm, ...
          'InitialMagnification', 400);

%  q
  % xlabel('x (mm)')
  %ylabel('depth (mm)')
%  axis on

  
% alphaIm = 0.8*(redIm > 0);
 %hold on
  %imagesc(xAxSel_mm, depthAxSel_mm, redIm(:,:,q), 'alphaData',
  %alphaIm); 
  %magesc(redIm, 'alphaData', alphaIm); 
  %rawnow
  %old off
  
  if 0
  figure(1)
  clf
  imagesc(showImScale)
  figure(2)
  clf
  imagesc()
  colorbar
  end
  
  drawnow
%  pausede
end

dt_s = min(diff(out.t_s));
ti_s = out.t_s(1):dt_s:out.t_s(end);

flowSumInterp = interp1(out.t_s, out.flowSum, ti_s, 'pchip', 0);

fsSig_Hz = 1/dt_s;
tAx_s = ti_s;
tuneParmsUS.approxPeriodFrac = 0.6; % force separation of peaks by
                                    % this fraction of period
tuneParmsUS.risingFitPolyOrder = 3;
tuneParmsUS.footGradThresh = 0.25;
tuneParmsUS.periodWindowFracReverse = 3/4;
tuneParmsUS.periodWindowFracForward = 1/3;
tuneParmsUS.risingFitGradWindow = 15e-3;
tuneParmsUS.riseTimeFrac = [0.05 0.95];
tuneParmsUS.maxMaxGradSpacing = 25e-3; % max seconds between max
                                       % grad
tuneParmsUS.lpf.periodScale = 0.5;
tuneParmsUS.lpf.n = 3;
tuneParmsUS.pressureLowLim = -Inf;
tuneParmsUS.pressureHighLim = Inf;
tuneParmsUS.fCutHigh = 0.5; % Hz
tuneParmsUS.fCutLowPeriod_Hz = 3; % Hz
tuneParmsUS.minFallTimeFactor = 1/20; % s
tuneParmsUS.tooCloseToPeakFracPeriod = 1/20;
%tuneParmsUS.fWall_Hz = 80;
tuneParmsUS.peakWinLenInd = 50; % how much of peak to leave behind
tuneParmsUS.periodTol = 0.25; 
tuneParmsUS.usPeakPrc = 66;
tuneParmsUS.usPeakThresh = 0.60; % to detect peak as frac of thresh
sigPre = flowSumInterp;
sigHPF = sigPre-mean(sigPre);
fCutLowPeriodNorm = tuneParmsUS.fCutLowPeriod_Hz/(fsSig_Hz/2); % Hz
[BL,AL] = butter(4, fCutLowPeriodNorm, 'low');
sigBPF = filtfilt(BL,AL, sigHPF);
approxPeriod = estimateperiod(sigBPF, fsSig_Hz);

tuneParmsUS.minRiseTime = approxPeriod*tuneParmsUS.minFallTimeFactor;
tuneParmsUS.lpf.fCut = 1/((approxPeriod*tuneParmsUS.lpf.periodScale)*...
                                fsSig_Hz/2);
tuneParmsUS.tooCloseToPeakTime = approxPeriod*...
    tuneParmsUS.tooCloseToPeakFracPeriod;  

sigGrad = gradient(sigBPF);
sig2Grad = gradient(sigGrad);
      
usPeakVal = prctile(sigGrad, tuneParmsUS.usPeakPrc);
peakGradThresh = tuneParmsUS.usPeakThresh*usPeakVal;
usPeaksInd0 = find( (sigGrad > peakGradThresh) & ...
                    sig2Grad > 0);

% repeat until we have all periodic peaks

usPeaksInd = usPeaksInd0;
approxPeriodFrac = tuneParmsUS.approxPeriodFrac;
approxPeriodInd = time2ind(approxPeriod,fsSig_Hz);
numDiscretePeaksOld = 0;

sigGrad0 = gradient(sigBPF);

footMark = [];
if isempty(footMark)
  [risingTimes, risingInd, sigGrad] = ...
      findperiodicpeaksiterfn(sigGrad0, peakGradThresh, ...
                              fsSig_Hz, ...
                              approxPeriod, tuneParmsUS);
        
  sigGradR0 = fliplr(sigGrad);
  [risingTimesBR, risingIndBR, sigGradBR] = ...
      findperiodicpeaksiterfn(sigGradR0, peakGradThresh, ...
                              fsSig_Hz, ...
                              approxPeriod, tuneParmsUS);
  
  risingIndR = length(sigGrad)-risingIndBR;
  risingTimesR = tAx_s(risingIndR);
  sigGradR = fliplr(sigGradBR);
  
  risingIndComb = [risingInd risingIndR];
  risingIndComb = unique(risingIndComb);
  
  footMark = tAx_s(risingIndComb);
  
  figure(3)
  clf
  subplot(3,1,1)
  plot(tAx_s, sigGrad)
  hold on
  plot(tAx_s(usPeaksInd), sigGrad(usPeaksInd),'.')
  plot(tAx_s(risingInd), sigGrad(risingInd),'g*')
  line(xlim, peakGradThresh* [1 1]);
  title('gradient');
  hold off
  subplot(3,1,2)
  plot(tAx_s, sigGradR)
  hold on
  %plot(tAx_s(usPeaksInd), sigGrad(usPeaksInd),'.')
  plot(tAx_s(risingIndR), sigGradR(risingIndR),'g*')
  line(xlim, peakGradThresh* [1 1]);
  title('gradient R');
  hold off
  subplot(3,1,3)
  plot(tAx_s, sigGradR)
  hold on
  %plot(tAx_s(usPeaksInd), sigGrad(usPeaksInd),'.')
  plot(tAx_s(risingIndComb), sigGradR(risingIndComb),'g*')
  line(xlim, peakGradThresh* [1 1]);
  hold off
  
  %      end
end

numMarks = length(footMark);
peakInd = [];
for r = 1:numMarks-1
  [mx, mxInd] = max(sigBPF(risingIndComb(r):risingIndComb(r+1)));
  peakInd(r) = mxInd + risingIndComb(r) - 1;
end

figure(4)
clf
plot(tAx_s, sigPre, 'r');
hold on
for r = 1:numMarks-1
  line(tAx_s(peakInd(r))*[1 1], ylim);
end

out.maxFlowImage = [];
out.maxFlowBImage = [];
out.peakTimeInterp_s = tAx_s(peakInd);

for r = 1:numMarks-1
  [dum, ind] = findclosestinvec(out.t_s, out.peakTimeInterp_s(r));
  out.maxFlowImage(:,:,r) =  out.flowImSto(:,:, ind);
  out.maxFlowBImage(:,:,r) =  out.BImSto(:,:, ind);  
end


segSet = [];
segImage = out.maxFlowImage;
for q =  1:numMarks-1
  
  % remove other flow sources
  cmPre_x = sum(BImStruct.xAxInterp_mm.*sum(segImage(:,:,q),1))/...
         sumall(segImage(:,:,q));
  cmPre_z = sum(BImStruct.zAxInterp_mm.*sum(segImage(:,:,q),2).')/...
         sumall(segImage(:,:,q));
  
  R2 = (BImStruct.XI_mm-cmPre_x).^2 + (BImStruct.ZI_mm-cmPre_z).^2;
  
  ind = find(R2 > tuneParms.maxDia_mm^2);

  thisFlow = segImage(:,:,q);
  thisFlow(ind)=0;
  
  cm_x_mm = sum(BImStruct.xAxInterp_mm.*sum(thisFlow,1))/...
         sumall(thisFlow);
  cm_z_mm = sum(BImStruct.zAxInterp_mm.*sum(thisFlow,2).')/...
         sumall(thisFlow);
    
  xz_mm = [cm_x_mm cm_z_mm];
  tuneParms.r0_mm = xz_mm(:);
  tuneParms.r0Min_mm = tuneParms.r0_mm-tuneParms.maxDia_mm/...
      tuneParms.r0DeltaAsFracMaxDia;
  tuneParms.r0Max_mm = tuneParms.r0_mm+tuneParms.maxDia_mm/...
      tuneParms.r0DeltaAsFracMaxDia;
  imParms = [];
  imParms.xAxInterp_mm = BImStruct.xAxInterp_mm;
  imParms.zAxInterp_mm = BImStruct.zAxInterp_mm;
  imParms.BIm = thisFlow;

  [prmOpt,  prmMax, prmMin, ellOptIm] = fitellmaskartertyfn(xz_mm, imParms, ...
                                                  tuneParms, plt);

  [ellIm] = prmOpt(5)*drawellipsefn(prmOpt(1:2), prmOpt(3), prmOpt(4), ...
                          BImStruct.XI_mm, BImStruct.ZI_mm);

  duplexIm = out.maxFlowBImage(:,:,q) + ellIm/2 + thisFlow/2; %out.flowImSto(:,:,q);
  figure(1)
  clf
  imshow( duplexIm, [0 256], 'colorMap', colMap, 'xdata', ...
         BImStruct.xAxInterp_mm, 'ydata', BImStruct.zAxInterp_mm,...
          'InitialMagnification', 400);

  duplexImFlow = out.maxFlowBImage(:,:,q) + thisFlow; %out.flowImSto(:,:,q);
  figure(2)
  clf
  imagesc(BImStruct.xAxInterp_mm, BImStruct.zAxInterp_mm, duplexImFlow);
  drawnow
  axis equal

  segSet.prmOpt(:,q) = prmOpt;
end

aDia_mm =  segSet.prmOpt(3,:)*2;
bDia_mm =  segSet.prmOpt(4,:)*2;

tuneParms.diaPercLimit = 50;

indSel = find(aDia_mm >= prctile(aDia_mm, tuneParms.diaPercLimit) & ...
           bDia_mm >= prctile(bDia_mm, tuneParms.diaPercLimit));
segSet.BImStruct = BImStruct;
segSet.prmOptSel = segSet.prmOpt(:,indSel);
segSet.maxFlowBImage = out.maxFlowBImage;
segSet.maxFlowImage = out.maxFlowImage
segSet.tuneParms = tuneParms;
segSet.imParms = imParms;
segSet.prmMax = prmMax;
segSet.prmMin = prmMin;
segSet.srcScript = srcScript;
segSet.h5File = h5File ;
segSet.h5FileCDI = h5FileCDI;
segSet.fifoLen = fifoLen;
segSet.d_mm = d_mm;

outFile = fullfile(pth, [mfile '_' fln '.mat']);
save(outFile, 'segSet');
lslrt(outFile)