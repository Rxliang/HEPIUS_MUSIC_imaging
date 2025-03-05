% seems deprecated. tried to read iq. think readmmodever1 is preferable.


vtHome = getenv('MUSIC_HOME');
vtData = getenv('MUSIC_DATA');

dataPath = fullfile(vtData, 'procwtfn');

overWrite = 1;

frameStart = 1;

pltB=0;
pltM=1;

% this is freq used to set PData vals
baseFreq_MHz = 2.8409;
baseLambda_mm = 1540/baseFreq_MHz/1e3;

if 0
  inFilePre = 'procwtfn_20190927_1702'
  fileNo = 4; 
  
  xLumenOffset_mm = 0; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 1;
  
  frameStart = 1;
  plt=0;
end

if 0
  inFilePre = 'procwtfn_20190927_2325'
  fileNo = 5; 
  
  xLumenOffset_mm = 0; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 3;
  zHalfSpan_mm = 1;
  
  frameStart = 1;
  plt=0;
end

if 0
  inFilePre = 'procwtfn_20190927_2357'
  fileNo = 5; 
  fileNo = 4; 
  fileNo = 3; 
  
  
  xLumenOffset_mm = -0.5; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 3;
  zHalfSpan_mm = 2;
  
  frameStart = 1;
  plt=0;
end


if 0

  inFilePre = 'procwtfn_20190928_0007'
  % to left of Doppler patch, 2 to 5
  fileNo = 2; 
  % to right of Doppler patch
  fileNo = 5; 
  
  % purposely off flow
  fileNo = 14; 
  
%  fileNo = 19; % higher pressure 
%  fileNo = 18; 
 % fileNo = 17; 
 % fileNo = 16;   
fileNo=14;
  xLumenOffset_mm = 0; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 2.5;
  zHalfSpan_mm = 2.5;
  
  frameStart = 1;
  plt=0;
end


if 0

  inFilePre = 'procwtfn_20191001_1101'
  fileNo=7;
  
  inFilePre = 'procwtfn_20191001_1119'
  fileNo=4;
  
  inFilePre = 'procwtfn_20191001_1220'
  fileNo=7;
  
  % good set
  fileNo=23; % var pressure 45mm
  xLumenOffset_mm = -1.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1.5;
  zHalfSpan_mm = 0;  
  
  fileNo=21; % var pressure 45mm
  xLumenOffset_mm = -1.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1.5;
  zHalfSpan_mm = 0;
  
  
  frameStart = 1;
  plt=0;
end

if 0

  inFilePre = 'procwtfn_20191001_1420'
  fileNo=7;
  
  inFilePre = 'procwtfn_20191001_1448'
  fileNo=17;
  
  inFilePre = 'procwtfn_20191001_1527'
  fileNo=6;
  
  inFilePre = 'procwtfn_20191001_2253'
  fileNo=14;
  
 inFilePre = 'procwtfn_20191001_2331'
  fileNo=11;
  
  inFilePre = 'procwtfn_20191002_0002'
  fileNo=2;
  
    inFilePre = 'procwtfn_20191002_0002'
  fileNo=10;

    inFilePre = 'procwtfn_20191002_0017'
  fileNo=4;
    fileNo=14;
    
    
    inFilePre = 'procwtfn_20191002_1416'
    fileNo=10;
    
    inFilePre = 'procwtfn_20191003_0017'
    fileNo=2;
  MRowAx = 'z';
  %MRowAx = 'x';
  
  
  xLumenOffset_mm = 2; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = -1.5; % offset from center of IQ matrix
  xHalfSpan_mm = 0.5;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
end
%inFile = [inFilePre '_' num2str(fileNo)];

if 0
  MRowAx = 'z';
%  MRowAx = 'x';
  
  xLumenOffset_mm = 0; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
%dataPath= fullfile(vtData, 'nicp', 'subject1_20191004183134/');

  sn = 'init_wt_s1_t1_sg15_i1';
  dataPath= fullfile(vtData, 'nicp', 'semaphoretest');

  inFile = sn;

end

if 0
  MRowAx = 'z';
%  MRowAx = 'x';
  
  xLumenOffset_mm = 0; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
  sn = 'subject1_wt_s4_t1_sg1_i1';
  
  subject1_20191008134056\subject1_wt_s4_t1_sg1_i1.im
  
  dataPath= fullfile(vtData, 'nicp', 'subject1_20191008112552');
  inFile = sn;

  
end


if 0
  MRowAx = 'z';
%  MRowAx = 'x';
  
  xLumenOffset_mm = 0; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
  sn = 'subject1_wt_s2_t1_sg1_i1';
  % retina at 32mm, applicator, 5,25,45 at 30s increments  
  dataPath= fullfile(vtData, 'nicp', 'subject1_20191008134056');
  inFile = sn;

  
end



if 0
  MRowAx = 'z';
%  MRowAx = 'x';
  
  xLumenOffset_mm = 0; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
  sn = 'subject1_wt_s3_t1_sg1_i1';
  % retina at 32mm, applicator, 5,25,45 at 30s increments  
  dataPath= fullfile(vtData, 'nicp', 'subject1_20191008133558');
  inFile = sn;

  
end

if 0
  MRowAx = 'z';
%  MRowAx = 'x';
  
  xLumenOffset_mm = 0; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
  sn = 'subject1_wt_s4_t1_sg1_i1';
  % retina at 32mm, applicator, 5,25,45 at 30s increments  
  dataPath= fullfile(vtData, 'nicp', 'subject1_20191008140349');
  inFile = sn;
  
end


if 1     
    % this is freq used to set PData vals
    baseFreq_MHz = 10.4167;
    baseLambda_mm = 1540/baseFreq_MHz/1e3;
    MRowAx = 'z';
%  MRowAx = 'x';
  
  xLumenOffset_mm = 0; %5.75; % offset from center of IQ matrix
  zLumenOffset_mm = 0; % offset from center of IQ matrix
  xHalfSpan_mm = 1;
  zHalfSpan_mm = 0.5;
  
  frameStart = 1;
  plt=0;
  sn = '20231120_2056_bwt_s1_t1_sg1_i1.rf';
  % retina at 32mm, applicator, 5,25,45 at 30s increments  
  dataPath= fullfile(vtData, 'nicp', 'sempahoretest');
                     
  inFile = sn;
  
end

startTime_s=0;
fileVersion=2;

if 1 | ~exist('iqSet', 'var') | isempty(iqSet)
  [outFile, iqSet] = readiqv2fn(inFile, startTime_s, dataPath, ...
                              fileVersion, overWrite);

  if isempty(iqSet)
    load(outFile);
  end
end

o_mm = iqSet.head.origin_wvl* baseLambda_mm;
xd_mm = iqSet.head.PDelta_wvl(1)* baseLambda_mm;
zd_mm = iqSet.head.PDelta_wvl(3)* baseLambda_mm;
sz = size(iqSet.IQMat);

xAx_mm = o_mm(1):xd_mm:o_mm(1)+(sz(2)-1)*xd_mm;
zAx_mm = o_mm(3):zd_mm:o_mm(3)+(sz(1)-1)*zd_mm;

if pltB
  figure(1)
  clf
  xlabel('x (mm)')
  zlabel('z (mm)')
end

xC_mm = mean(xAx_mm);
xLumen_mm = xC_mm + xLumenOffset_mm;
xL_mm = xLumen_mm - xHalfSpan_mm;
xH_mm = xLumen_mm + xHalfSpan_mm;

xIndM = find(xAx_mm <= xH_mm & xAx_mm >= xL_mm);

xIndML = find(xAx_mm <= xLumen_mm & xAx_mm >= xL_mm);
xIndMH = find(xAx_mm >= xLumen_mm & xAx_mm <= xH_mm);

[dum, xIndC] = findclosestinvec(xAx_mm, xLumen_mm);

zC_mm = mean(zAx_mm);
zLumen_mm = zC_mm + zLumenOffset_mm;
zL_mm = zLumen_mm - zHalfSpan_mm;
zH_mm = zLumen_mm + zHalfSpan_mm;

zIndML = find(zAx_mm <= zLumen_mm & zAx_mm >= zL_mm);
zIndMH = find(zAx_mm >= zLumen_mm & zAx_mm <= zH_mm);
zIndM = find(zAx_mm <= zH_mm & zAx_mm >= zL_mm);

zC_mm = mean([zL_mm zH_mm]);

if isempty(zIndM)
  [dum, zIndM] = findclosestinvec(zAx_mm, zLumen_mm);
end

%
%[dum, zIndC] = findclosestinvec(zAx_mm, zLumen_mm);

%MRowAx = 'z';
  %MRowAx = 'x';%zAxM_mm = zAx_mm(zIndC);

%szM = [length(xIndM)+1 sz(3)-frameStart+1]; % accommodate fwd lag one

if MRowAx == 'x'
  szM = [sz(2) sz(3)-frameStart+1]; % accommodate fwd lag
                                            % one
else
  szM = [sz(1) sz(3)-frameStart+1]; % accommodate fwd lag  
end

if 1 % ~exist('MMat', 'var') | isempty(MMat)
  MMat = zeros(szM, 'single');

  cnt = 1;
  for i = frameStart:size(iqSet.IQMat,3)
    thisB = iqSet.IQMat(:,:,i);
    if MRowAx == 'x'
      MMat(:,cnt) = mean(thisB(zIndM, :),1); %[xIndM
                                             %xIndM(end)+1]); 
    else
      MMat(:,cnt) = mean(thisB(:, xIndM), 2); 
    end
    
    cnt=cnt+1;
    
    if pltB
      imagesc(xAx_mm, zAx_mm, sqrt(abs(thisB)));
      xlabel('mm');
      ylabel('mm');
      if 1
        line(xLumen_mm*[1 1], ylim);
        hold on
        line(xlim, zLumen_mm*[1 1]);
        %      line(xlim, zL_mm*[1 1]);      
        %      line(xlim, zH_mm*[1 1]);            
        if MRowAx == 'x'
          line(xL_mm*[1 1], ylim);
          line(xH_mm*[1 1], ylim);      
        else
          line(zL_mm*[1 1], ylim);
          line(zH_mm*[1 1], ylim);      
        end
        
        
        
      end
      drawnow
      %    pause(0.01)
    end
  end
  
end

if MRowAx == 'x'
   MMatSel = MMat([xIndML xIndMH],:);
   rowAx_mm = xAx_mm;
else
   MMatSel = MMat([zIndML zIndMH],:).';
   rowAx_mm = zAx_mm;
end

%MMatSel = MMat([xIndML],:);

peakInd = [];
span = 0;
refInd = round(size(MMatSel,2)/2);


for i = 1:(size(MMatSel,2)- span)
  m1 = abs(MMatSel(:,refInd));
  m1 = m1-mean(m1);
  
  m2 = abs(MMatSel(:,i+span));
  m2 = m2-mean(m2);
  
  [r, lags] = xcorr(m2,m1);
  r0Ind = find(lags==0);
  [dum, peakIndPre] = max(r);
  
  peakInd(i) = peakIndPre-r0Ind;
  
  if 0
    figure(10)
    plot(r)
    pausede
  end
  
end

if MRowAx == 'x'
  diff_mm = (xAx_mm(2)-xAx_mm(1)); 
else
  diff_mm = (zAx_mm(2)-zAx_mm(1)); 
end


d2s = 24*3600; 
tAx_s = d2s*iqSet.head.datenum;
tAx_s = tAx_s-tAx_s(1);

figure(10)
clf
plot(tAx_s(1:length(peakInd)), -peakInd*diff_mm)

if MRowAx == 'x'
  phiL = c3mfn(MMat.', xIndML);
  phiH = c3mfn(MMat.', xIndMH);
else
  phiL = c3mfn(MMat.', zIndML);
  phiH = c3mfn(MMat.', zIndMH);
end
  
sL_mm = cumsum(phiL);
sH_mm = cumsum(phiH);


figure(2)
clf
subplot(3,1,1)
plot(tAx_s(1:end-1), sL_mm)
subplot(3,1,2)
plot(tAx_s(1:end-1), sH_mm)
subplot(3,1,3)
plot(tAx_s(1:end-1),-sH_mm+sL_mm)


if pltM
  figure(3)
  clf
  %freq_MHz=4.6;
  
  freq_MHz= iqSet.head.freq_MHz;
  
%  alpha=0.54;
  alpha=1.8;
  
  imM = abs(MMat);
  
  if MRowAx == 'z'
    [atten, comp] = usattenfn(rowAx_mm, freq_MHz, alpha);
    imM = imM .* repmat(comp(:),1, size(MMat,2));
  end
 
  imagesc(tAx_s, rowAx_mm, imM)
  hold on
  if MRowAx == 'x'
    line(xlim, xLumen_mm*[1 1], 'color', 'r', 'linestyle', '-.');
    line(xlim, xL_mm*[1 1], 'color', 'r', 'linestyle', '-.');
    line(xlim, xH_mm*[1 1], 'color', 'r', 'linestyle', '-.');  
    ylabel('x (mm)')
  else
    line(xlim, zLumen_mm*[1 1], 'color', 'r', 'linestyle', '-.');
    line(xlim, zL_mm*[1 1], 'color', 'r', 'linestyle', '-.');
    line(xlim, zH_mm*[1 1], 'color', 'r', 'linestyle', '-.');  
    ylabel('z (mm)')
  end
    
  xlabel('t (s)')

  hold off
end

outFile = fullfile(dataPath, [sn '.mat']);
save(outFile, 'tAx_s', 'rowAx_mm', 'imM');
lslrt(outFile)


