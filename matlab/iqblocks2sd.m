mfile = mfilename;
vtData = getenv('MUSIC_DATA');

inPath = fullfile(vtData, '0627MUSICPIG');
matFileIn = 'FUS_pig_afterkataminepreFUS20240627_194601_sdl_s1_sg1.mat';
matFile = fullfile(inPath, matFileIn);
lslrt(matFile)
s = load(matFile, 'sdOut');

Nwind = s.sdOut.SDop.nIqfft;

w = hammingvs(Nwind);
w = w/sum(w);
%specWindow = w(:)*ones(1,s.sdOut.sgp.R);
        
% multiply all cols of IQ by window and apply DFT to each col
sfBlock = abs(fft(s.sdOut.IQSlideBlocks.*w, s.sdOut.SDop.nfft)).^2;

if isequal(s.sdOut.SDop.revFlowDir, 1)
  sfBlock = flipud(sfBlock);
end

prcOutlierCutPerc = 99;
sfBlockFix = sfBlock;
prcThresh = prctile(sfBlockFix(:), prcOutlierCutPerc);
indKill = find(sfBlockFix > prcThresh);
sfBlockFix(indKill) = NaN;

% these are not even
iiCell = subframeInd100fn(s.sdOut.sgp, s.sdOut.SDop.PWdopPRF)


Ts = 

tAx = 

figure(2)
clf
imagesc(sfBlockFix)
axis xy




  indKill = find(medRows > prcThresh);
  sfBlockFix(:, indKill) = NaN;


% outlier removal
indKillOld = [];
while 1
  medRows = median(sfBlockFix,1);
  indKill = find(medRows > prcThresh);
  sfBlockFix(:, indKill) = NaN;
  
  if isempty(indKill) %| length(indKill) == length(indKillOld)
    break
  end
  
%  indKillOld = indKill;
%  indKillSelInd = find(indKill > 1);
%  theseCols= indKill(indKillSelInd)
%  sfBlockFix(:, theseCols) = sfBlockFix(:, theseCols-1);
end


%sfBlockFix(:, indKill) = NaN;
medVal = median(sfBlockFix,1);




% was 30, but surely should not show the pulse level shift
NumNoiseHist = 300; %noise floor estimation median window size

%normalizing spec display data:
%calc noise floor
% sfBlock has R columns, so this looks at last column

%noiseFloorHistory = [median(sfBlock(:,end),1),noiseFloorHistory];
%noiseFloorHistory(NumNoiseHist:end)=[];

medVal = median(sfBlock,1);

% try to use a similar number of samples for medial filtering
medFiltWinLen = NumNoiseHist * s.sdOut.sgp.R;

medValFilt = medfilt2(medVal, [1 medFiltWinLen], 'symmetric');

    nfhMid = median(noiseFloorHistory);
    if isempty(noiseFloor)
        noiseFloor = nfhMid;
    else
        noiseFloor = noiseFloor*SDop.noisePersist + ...
            nfhMid*(1-SDop.noisePersist);
    end

    sfBlock = sfBlock.'/noiseFloor; %normalize to noise floor

