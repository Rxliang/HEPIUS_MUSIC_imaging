
vtHome = getenv('BNS_HOME');
vtData = getenv('BNS_DATA');

NF = size(imSet.im,3);
figure(1)
clf
%indSel=120:500;
%indSel= 93:500 % vein

tAxPre_s = 3600*24*(imSet.datenum- min(imSet.datenum));
[mn, indMin] = min(tAxPre_s);
[mx, indMax] = max(tAxPre_s);

indSelPre = [indMin:length(tAxPre_s)]; % 1:indMax];

tAx_s = tAxPre_s(indSelPre);
dt_s = diff(tAx_s);

minDiff = min(dt_s);
maxDiff = max(dt_s);
midDiff = mean([minDiff maxDiff]);

indAx = 1:length(dt_s);
ind = 1+find((dt_s > midDiff) & indAx >= P.numRepeats-1, 1);

colSel = 1:633; %size(imSet.im,2); %518;
colSel = 1:size(imSet.im,2); %518;
%colSel = 200:550;

indSel = indSelPre(ind:end);

imSelPre = imSet.im(:,colSel, indSel);

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

colormap(gray)
for q = 1:length(indSel)
  imagesc(imSelPre(:,:,q).^0.5);
  q
  drawnow
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

figure(2); clf;
imDiffFilt = [];
for q = 1:size(imDiff,3)
  %imagesc(imDiff(:,:,q))
  imDiffFilt(:,:,q) = imgaussfilt(imDiffMean(:,:,q),3);
  imagesc(imDiffFilt(:,:,q))
  drawnow
  
end

pwrAllImDiff = squeeze(mean(imDiffFilt,[1 2]));
prctTime = prctile(pwrAllImDiff, 85);
highFlowInd = find(pwrAllImDiff >= prctTime);
highFlowImDiffMean = mean(imDiffMean(:,:,highFlowInd),3);
prctSpace = prctile(highFlowImDiffMean(:), 95.5);


%for q = 1:size(imDiff,3)
%  pwrImDiff = meanall(imDiffMean(:,:,q));
%end
  
imBMean = mean(imB,4);
colMap = grayscaleCPAmap;

for q = 1:size(imDiff,3)
  ind = find(imDiffFilt(:,:,q) >= prctSpace);
  thisFlowIm = imDiffFilt(:,:,q);
  
%  indnan = find(isnan(thisFlowIm));
%  thisFlowIm(indnan)= 0;
  thisBIm = imBMean(:,:,q).^0.5;  
  mnFlow = meanall(thisFlowIm);
  mnB = meanall(thisBIm);
  flowImScaled = thisFlowIm*mnB/mnFlow;
  showIm = thisBIm;
  showImScale = (showIm-minall(showIm))/(maxall(showIm)-minall(showIm))*127;
  %imagesc(showIm, [minall(thisBIm) maxall(thisBIm)]);
  imshow(showImScale, [0 256], 'colorMap', colMap)%
 % imshow(showImScale%
  
  %rgbIm = repmat(zeros(size(showIm)),1,1,3);
  redIm = zeros(size(showIm));
  fmn = minall(flowImScaled(ind));
  fmx = maxall(flowImScaled(ind));
  redIm(ind) = (flowImScaled(ind)-fmn)*128/fmx + 128;
  %/prctSpace;
  %rgbIm(:,:,1) =  redIm;
  alphaIm = 0.8*(redIm > 0);
%  if ~isempty(ind)
 %   return
  %end
  hold on
  %image(rgbIm, colMap, 'alphaData', alphaIm);
  imagesc(redIm, 'alphaData', alphaIm); 
%  imshow(rgbIm/maxall(showIm)) %, [minall(thisBIm) maxall(thisBIm)]);  
  q
  drawnow
  hold off
  pause(0.05)
end





















if 0
outFile = fullfile(vtData, [mfile '_' num2str(imSet.datenum(1)) ...
                    '.mat']);

save(outFile, 'P', 'imSet');
lslrt(outFile)
end
