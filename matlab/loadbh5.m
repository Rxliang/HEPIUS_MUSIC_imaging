mfile = mfilename;

vtData = getenv('MUSIC_DATA');
srcScript = 'bflowmusicl818';
h5File = fullfile(vtData, srcScript, 'bflowmusicl818_20231202_1815.h5'); % good radial
h5File = fullfile(vtData, srcScript, 'bflowmusicl818_20231202_1827.h5');
h5File = fullfile(vtData, srcScript, 'bflowmusicl818_20231202_1835.h5');

[pth, fln, ext] = fileparts(h5File);

BFrame = h5read(h5File, '/data/B');
datenumFrame = h5read(h5File, '/data/datenumFrame');

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

figure(1)
clf
for q = 1:size(BFrame,3)
  
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

  pausede  
end
