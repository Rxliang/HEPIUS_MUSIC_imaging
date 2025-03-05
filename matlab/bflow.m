mfile = 'pseudobflow'
vtHome = getenv('MUSIC_HOME');
vtData = getenv('MUSIC_DATA');

if 0
  load(fullfile(vtData,'gel18rylns256bf_737845.0444.mat'));
  site = 'carotid';
end

if 0
  load(fullfile(vtData,'gel18rylns256bf_737845.0366.mat'))
  site = 'carotid_10MHz'
end

if 1
  load(fullfile(vtData,'bflow2020/gel18rylns256bf_737844.9986.mat'))
  site = 'radial_distal_to_bifurcation_10MHz'
end

makeMovie = 1;

if 1
  movieSpec.movieName = [vtData '/' mfile '_' site];
  vidObj = VideoWriter(movieSpec.movieName);
  vidObj.FrameRate = 15;
  open(vidObj);
end

lambda_mm = 1e-3*1540/Trans.frequency;

parms.BMaskPerctile = 75;

if 1 % carotid
  parms.BMaskPerctile = 75;
  parms.BMaskPerctile = 90;  
  parms.frameStart = 1;
  parms.frameStart = 1;
  flowTimePercentile = 85;
  flowSpacePercentile = 87.5;  
end

if 0 % radial
  parms.BMaskPerctile = 75;
  parms.frameStart = 1;
  parms.frameStart = 1;
  flowTimePercentile = 90;
end
NF = size(imSet.im,3);
figure(1)
clf
%indSel=120:500;
%indSel= 93:500 % vein

depthAx_wvl = PData.Origin(3):PData.PDelta(3):PData.Origin(3)+(PData.PDelta(3)*(PData.Size(1)-1));
depthAx_mm = depthAx_wvl*lambda_mm;

xAx_wvl = PData.Origin(1):PData.PDelta(1):PData.Origin(1)+(PData.PDelta(1)*(PData.Size(2)-1));
xAx_mm = xAx_wvl*lambda_mm;

[X_mm,Z_mm]= meshgrid(xAx_mm, depthAx_mm);

tAxPre_s = 3600*24*(imSet.datenum- min(imSet.datenum));
[mn, indMin] = min(tAxPre_s);
[mx, indMax] = max(tAxPre_s);

indSelPrePre = [indMin:length(tAxPre_s)]; % 1:indMax];
indSelPre =  indSelPrePre(parms.frameStart:end);

tAx_s = tAxPre_s(indSelPre);
dt_s = diff(tAx_s);

minDiff = min(dt_s);
maxDiff = max(dt_s);
midDiff = mean([minDiff maxDiff]);

indAx = 1:length(dt_s);
ind = 1+find((dt_s > midDiff) & indAx >= P.numRepeats-1, 1);


colSel = 1:633; %size(imSet.im,2); %518;
colSel = 1:size(imSet.im,2); %518;
colSel = 1:size(imSet.im,2); %518;
colSel = 100:700;

rowSel = 1:size(imSet.im,1); %518;

if strcmp(site(1:4), 'carotid')
  rowSel = 1:650;
end

%colSel = 300:size(imSet.im,2); %518;

indSel = indSelPre(ind:end);

XSel_mm = X_mm(rowSel, colSel);
ZSel_mm = Z_mm(rowSel, colSel);
xAxSel_mm = xAx_mm(colSel);
depthAxSel_mm = depthAx_mm(rowSel);

imSelPre = imSet.im(rowSel,colSel, indSel);

if 0
  % need to align ensemble
  ensemInd = 1:P.numRepeats:NF;
  indEnsemStart = find(ensemInd >= indMin);
  indSelStart = find(indSelPre>=ensemInd(indEnsemStart(1)))
  indSel = indSelPre(indSelStart:end);
  
  indSelPre = 171:NF;
  indSelPre = 103:NF;
  indSelPre = 90:NF;
  indSelPre = 92:NF;
  
  % need to align ensemble
  ensemInd = 1:P.numRepeats:NF;
  indEnsemStart = find(ensemInd >= indSelPre(1));
  indSelStart = find(indSelPre>=ensemInd(indEnsemStart(1)))
  indSel = indSelPre(indSelStart:end);
end

if 0
colormap(gray)
for q = 1:length(indSel)
  imagesc(imSelPre(:,:,q).^0.5);
  q
  drawnow
end
end

lenSel = length(indSel);
indDiffPre = [];
lenDiff=[];
for p = 1:P.numRepeats
  indDiffPre{p} = p:P.numRepeats:lenSel;
  lenDiff(p) = length(indDiffPre{p});
end

indDiff =[];
for p = 1:P.numRepeats
  indDiff(p,:) = indDiffPre{p}(1:min(lenDiff));
end

diffRowOffset = 2:1:P.numRepeats;
numAbsDiffs = length(diffRowOffset);

imDiff = [];
imB = [];
for q = 1:numAbsDiffs
  imSubtrahend = imSelPre(:,:, indDiff(diffRowOffset(q),:));
  imMinuend = imSelPre(:,:, indDiff(diffRowOffset(q)-1,:));
  
  imDiff(:,:,:,q) =  abs(imMinuend -  imSubtrahend);
  imB(:,:,:,q) =  (imMinuend +  imSubtrahend)/2;  
    
end

imDiffMean = mean(imDiff,4);
imBMean = mean(imB,4);

figure(2); clf;
imDiffFilt = [];
for q = 1:size(imDiff,3)
  %imagesc(imDiff(:,:,q))

  thisBIm = imBMean(:,:,q).^0.5;  
  indMaskB = find(thisBIm > prctile(thisBIm(:), parms.BMaskPerctile) | ZSel_mm < 2 | ...
                  XSel_mm <  (XSel_mm(1,1)+1) | XSel_mm > (XSel_mm(1,end))-1);
  
  thisBMask = zeros(size(thisBIm));
  thisBMask(indMaskB) =1;
  
  imFilt = imgaussfilt(imDiffMean(:,:,q),7);
  imFilt(indMaskB)=0;
  imDiffFilt(:,:,q) = imFilt;
  
  imagesc(imDiffFilt(:,:,q))
  drawnow
end

pwrAllImDiff = squeeze(mean(imDiffFilt,[1 2]));
prctTime = prctile(pwrAllImDiff, flowTimePercentile);
highFlowInd = find(pwrAllImDiff >= prctTime);
highFlowImDiffMean = mean(imDiffMean(:,:,highFlowInd),3);
prctSpace = prctile(highFlowImDiffMean(:), flowSpacePercentile);

%for q = 1:size(imDiff,3)
%  pwrImDiff = meanall(imDiffMean(:,:,q));
%end
  

colMap = grayscaleCPAmap;

sto.redIm = [];
sto.showImScale = [];
for q = 1:size(imDiff,3)
  ind = find(imDiffFilt(:,:,q) >= prctSpace);
  thisFlowIm = imDiffFilt(:,:,q);
   thisBIm = imBMean(:,:,q).^0.5; 
%  indnan = find(isnan(thisFlowIm));
%  thisFlowIm(indnan)= 0;
  mnFlow = meanall(thisFlowIm);
  mnB = meanall(thisBIm);
  flowImScaled = thisFlowIm*mnB/mnFlow;
  flowImScaled(indMaskB)=0;
  showIm = thisBIm;
  showImScale = (showIm-minall(showIm))/(maxall(showIm)-minall(showIm))*127;
  %imagesc(showIm, [minall(thisBIm) maxall(thisBIm)]);
 % imshow(showImScale, [0 256], 'colorMap', colMap)%
 % imshow(showImScale%
  
  %rgbIm = repmat(zeros(size(showIm)),1,1,3);
  redIm = zeros(size(showIm));
  fmn = minall(flowImScaled(ind));
  fmx = maxall(flowImScaled(ind));
  redIm(ind) = (flowImScaled(ind)-fmn)*128/(fmx-fmn) + 128;
  %/prctSpace;
  %rgbIm(:,:,1) =  redIm;
%  alphaIm = 0.8*(redIm > 0);
%  if ~isempty(ind)
 %   return
  %end
 
%  imshow(rgbIm/maxall(showIm)) %, [minall(thisBIm) maxall(thisBIm)]);  
  q
%  drawnow
%  hold off
 % pause(0.025)
  sto.showImScale(:,:,q)=showImScale;
  sto.redIm(:,:,q) = redIm;
end

figure(1)
clf
for q = 1:size(imDiff,3)
  imshow(sto.showImScale(:,:,q), [0 256], 'colorMap', colMap, 'xdata', ...
         xAxSel_mm, 'ydata',depthAxSel_mm); 
  xlabel('x (mm)')
  ylabel('depth (mm)')
  axis on

  
  alphaIm = 0.8*(sto.redIm(:,:,q) > 0);
  hold on
  %image(rgbIm, colMap, 'alphaData', alphaIm);
  %imagesc(sto.redIm(:,:,q), 'alphaData', ...
   %       alphaIm); 
  imagesc(xAxSel_mm, depthAxSel_mm, sto.redIm(:,:,q), 'alphaData', alphaIm); 
  drawnow
  hold off
  if makeMovie
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end
  
  %pause(0.025)
end

if makeMovie
  close(vidObj);
end













if 0
outFile = fullfile(vtData, [mfile '_' num2str(imSet.datenum(1)) ...
                    '.mat']);

save(outFile, 'P', 'imSet', 'Trans', 'TW', 'PData', 'TPC', 'SeqControl', ...
     'parms', '-v7.3');
lslrt(outFile)
end
