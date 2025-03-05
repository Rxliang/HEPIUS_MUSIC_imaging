% read multiframe m-mode file from procrflinefn post 2019/10/21

vtHome = getenv('MUSIC_HOME');
vtData = getenv('MUSIC_DATA');

% added x,z for inst one to header, as well as freq
inFile = ['D:\vittamed_data\nicp\semaphoretest\' ...
          'init_bwt_s1_t1_sg1_i1_fix11.rf'];
inFile = 'D:\vittamed_data\nicp\subject1_20191026182221\subject1_bwt_s1_t1_sg1_i1.rf';

if 1
    %  sn = '20231120_2056_bwt_s1_t1_sg1_i1.rf';
    %      sn = '20231120_2131_bwt_s1_t1_sg2_i1.rf';
          sn = '20231120_2200_bwt_s1_t1_sg1_i1.rf'; %  high gain, high zin
                                                    %          sn = '20231120_2220_bwt_s1_t1_sg1_i1.rf';
          sn = '20231120_2220_bwt_s1_t1_sg1_i1.rf'; %  low gain, low zin,
          sn = '20231120_2233_bwt_s1_t1_sg1_i1.rf'; %  low gain, low zin, raise V in middle
          sn = '20231120_2248_bwt_s1_t1_sg1_i1.rf'; %  low gain, low zin, raise V in middle, zero tx apod for all seqs
                    sn = '20231120_2257_bwt_s1_t1_sg1_i1.rf'; %  low gain, low zin, raise V in middle, zero tx apod for all seqs
  % retina at 32mm, applicator, 5,25,45 at 30s increments  
  dataPath= fullfile(vtData, 'nicp', 'semaphoretest');
  inFile = fullfile(dataPath, sn);
end

lslrt(inFile)

MSpec = [];
MSpec.startDepthCorr_mm = 0;
MSpec.tissueAtten_dBpMHz = 1.0;
MSpec.tSel_s = [0 Inf];
MSpec.plt=10;
MSet = readrfbwtfn(inFile, MSpec);
title(titstrfn(['M-mode: ' sn]));
colormap(gray)

figFile = [inFile '.png'];
print('-dpng', figFile)
lslrt(figFile)








return


deltaDepth_mm = 0.75;

lumenDepth_mm = 53.57% z_mm;
tSel_s = [35 50];

lumenDepthWt_mm = 46.86% z_mm;
tSel_s = [35 50];

%lumenDepth_mm = 38.52% z_mm;
%tSel_s = [10 18];

deltaDepth_mm = 1.5;
lumenDepthWt_mm = 15.8% z_mm;
tSel_s = [50 70];

deltaDepth_mm = 1;
lumenDepthWt_mm = 50;% z_mm;
tSel_s = [20 40];


deltaDepth_mm = 2;
%lumenDepthWt_mm = 6 %zGate_mm; % z_mm;
tSel_s = [5 30];

distalDepth_mm = lumenDepthWt_mm-deltaDepth_mm;
proximalDepth_mm = lumenDepthWt_mm+deltaDepth_mm;

colSel = find(tAx_s >= tSel_s(1) & tAx_s <= tSel_s(2));
MMat = M(:, colSel);

depthAx_mm = zAx_mm;
distalInd = find(depthAx_mm <=  lumenDepthWt_mm & ...
                 depthAx_mm >=  distalDepth_mm);
proximalInd = find(depthAx_mm >=  lumenDepthWt_mm & ...
                   depthAx_mm <=  proximalDepth_mm);

line(xlim, lumenDepthWt_mm*[1 1], 'linestyle','--','color', 'y');
hold off

MPart = MMat([distalInd; proximalInd],:);

phiD = c3mfn(MMat.', distalInd);
phiP = c3mfn(MMat.', proximalInd);  
distensionD_mm = cumsum(phiD,2);
distensionP_mm = cumsum(phiP,2);  

  %distensionFilt_mm = detrend(distension_mm); %filtfilt(b, a,
  %distension_mm); 
distension_mm = +distensionD_mm+distensionP_mm;

tAxWT_s = tAx_s(colSel(1:end-1));
figure(11)
subplot(311)
hp = plot(tAxWT_s, detrend(distensionD_mm));
hp = plot(tAxWT_s, (distensionD_mm));
title('distal')

subplot(312)
%hp = plot(tAxWT_s, detrend(distensionD_mm));
hp = plot(tAxWT_s, (distensionP_mm));
title('proximal')

subplot(313)
%hp = plot(tAxWT_s, detrend(distensionD_mm));
hp = plot(tAxWT_s, (distension_mm));
%ylim([-0.5 0.5]);

%    filtD_mm = filtfilt(filtParms.BH,filtParms.AH, double(distensionD_mm));

 %   hp = plot(handleWT, tAxWT_s, filtD_mm);



    