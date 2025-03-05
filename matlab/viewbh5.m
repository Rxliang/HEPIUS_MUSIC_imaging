mfile = mfilename;


h5File = ['D:\vittamed_data\gel18rylns256bf\' ...
          'gel18rylns256bf_20200430_0047.h5'];

[pth, fln, ext] = fileparts(h5File);

if ~exist('segSet', 'var') | isempty(segSet)
  matFileIn = fullfile(pth, ...
                       'viewbh5_gel18rylns256bf_20200430_0047.mat');
  load(matFileIn, 'segSet');
  tuneParms = segSet.tuneParms;
  lambda_mm =  segSet.lambda_mm;
  prmOpt_mm = segSet.prmOpt_mm;
  raInd = segSet.raInd;
else
  
  raInd = [];
  mnInd = [];  
  raInd{1} = 534:669;
  mnInd{1} = 671:825;
  raInd{2} = 826:935;
  mnInd{2} = 979:1272;
  raInd{3} = 1277:1355;
  mnInd{3} = 1362:1459;

  lambda_mm = 0.1478;
  segSet.lambda_mm = lambda_mm;
  segSet.raInd = raInd;
  segSet.mnInd = mnInd;
  segSet.h5File = h5File;
  
  
  tuneParms.a_mm = 1.33;
  tuneParms.b_mm = 1.33;
  
  tuneParms.maxDia_mm = 5;
  tuneParms.aMax_mm = tuneParms.maxDia_mm*0.45;
  tuneParms.aMax_mm = 2;  
  tuneParms.bMax_mm = 2; % tuneParms.maxDia_mm*0.45;

  tuneParms.aMin_mm = 1.3;
  tuneParms.bMin_mm = 1.5;
  tuneParms.annuPrctile = 90;
  tuneParms.annuRatio = 4;
  tuneParms.lumenSampleHalfWin_mm = 0.5;
  tuneParms.maxDiaGrowFactor = 1.3;
  tuneParms.wallThickness_mm = 0.5;
end

outFile = fullfile(pth, [ mfile '_' fln '.mat']);

BFrame = h5read(h5File, '/data/B');

BPDelta = h5read(h5File, '/header/PDelta');
BOrigin = h5read(h5File, '/header/Origin');
BSize = h5read(h5File, '/header/Size');

xAxB_mm = (BOrigin(1):BPDelta(1):BOrigin(1)+BPDelta(1)*(BSize(2)-1))* ...
          lambda_mm;

zAxB_mm = (BOrigin(3):BPDelta(3):BOrigin(3)+BPDelta(3)*(BSize(1)-1))* ...
          lambda_mm;

dx_mm = diff(xAxB_mm(1:2));
dz_mm = diff(zAxB_mm(1:2));

d_mm = max([dx_mm dz_mm]);

xAxInterp_mm = xAxB_mm(1):d_mm:xAxB_mm(end);
zAxInterp_mm = zAxB_mm(1):d_mm:zAxB_mm(end);

[X_mm, Z_mm] = meshgrid(xAxB_mm, zAxB_mm);
[XI_mm, ZI_mm] = meshgrid(xAxInterp_mm, zAxInterp_mm);

if 0 
figure(1)
clf
for q = raInd{1}(2) % 540:size(BFrame,3)
  BIm = interp2(X_mm, Z_mm,  BFrame(:,:,q).^.35, XI_mm, ZI_mm);
  q
  imagesc(xAxInterp_mm, zAxInterp_mm, BIm)
  colormap(gray)
  axis equal
  drawnow

%  pausede
end
end

if ~isfield(segSet, 'badFitInd') | isempty(segSet.badFitInd)
  segSet.badFitInd = [];
end


figure(1)
clf
for q = raInd{2}(80:end) % 540:size(BFrame,3)
  %if length(prmOpt_mm) > q & ~isempty(find(prmOpt_mm(:,:,q)~=0))
  %  continue
  %end
  
  figure(1)
  clf
  BIm = interp2(X_mm, Z_mm,  BFrame(:,:,q).^.35, XI_mm, ZI_mm);
  q
  imagesc(xAxInterp_mm, zAxInterp_mm, BIm)
  colormap(gray)
  axis equal
  drawnow

  if 0 
  [xz_mm] = ginput(1);
  xz_mm = xz_mm.';
  end 
  
  % xz_mm = [6.6; 4.0];


  imParms.BIm = BFrame(:,:,q).^.55;
  imParms.xAxInterp_mm = xAxInterp_mm;
  imParms.zAxInterp_mm = zAxInterp_mm;

  while 1 
    
    if (size(prmOpt_mm,3) < q) | isempty(find(prmOpt_mm(:,:,q))~=0)
      tuneParms.r0_mm = [xz_mm(1); xz_mm(2)];
      tuneParms.r0Max_mm = tuneParms.r0_mm+[0.25; 0.25];
      tuneParms.r0Min_mm = tuneParms.r0_mm-[0.25; 0.25];
      figure(2)
      disp('Click center of artery, or space to reuse last')
      [x, z, but] = ginput(1);
%      x = 3.41; z =  3.91;
%      x = 4.75; z =  4.34;
      %xz_mm = [x z];

      but=1;
      if but == 1
        xz_mm = [x z];
      end
        
      plt.solve = 1;
    else
      plt.solve = 0;
      xz_mm = prmOpt_mm(1:2,:,q);
      plt.prmOpt_mm = prmOpt_mm(:,:,q); 
    end
    
    plt.flag = 1;
    plt.figNo = 2;
    
    if plt.solve
      disp('solving ...')
    end
    [prmOpt_mm(:,:,q), prmMax_mm, prmMin_mm, ellOpt] = fitellartertyfn(xz_mm, imParms, tuneParms, plt);
    [prmOpt_mm(:,:,q) prmMax_mm prmMin_mm];

    if plt.solve
      break
    end
    
    if ~plt.solve
      figure(2)
      disp('q - bad seg, space - ok, click to recalc')
      [x, z, but] = ginput(1);
    end
    
    badInd = find(segSet.badFitInd == q);
    if ~isempty(badInd)
      disp('bad fit already noted')
      %break
    end
    
    
    but
    if but == 'q'
      disp('bad fit noted')
      segSet.badFitInd = [segSet.badFitInd q];
      break
    end
    
    if but==32
      segSet.badFitInd = setdiff(segSet.badFitInd, q);
      break
    else
      xz_mm = [x z];
      prmOpt_mm(:,:,q) = prmOpt_mm(:,:,q)*0;
      disp('Will optimize on next loop')
    end
  end
  
end


segSet.tuneParms = tuneParms;
segSet.imParms = imParms;
segSet.prmOpt_mm = prmOpt_mm;
segSet.prmMax_mm = prmMax_mm;
segSet.prmMin_mm = prmMin_mm;

save(outFile, 'segSet');
lslrt(outFile)