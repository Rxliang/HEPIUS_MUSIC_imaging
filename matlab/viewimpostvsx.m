% view images after VSX runs  PlaneWaveBmode_PlaneWaveDoppler20231215

mfile = mfilename;
vtData = getenv('MUSIC_DATA');
dataPath = fullfile(vtData);

[dum, inMatFile, dum] = fileparts(inMat);

outFile = fullfile(dataPath, 'pig', [mfile '_' ...
                   inMatFile '.mat'])
imb = 3;

if imb ~=3
  pd = imb;
else
  pd = 1;
end

dmFreq_MHz = Receive(2).demodFrequency;
lambda_mm = 1e-3*P.c/dmFreq_MHz;
ax = makeimaxisfn(ImgData{imb}(:,:,1,1), lambda_mm*PData(pd).PDelta([1 3]));
ax.vAxB = ax.vAxB-ax.vAxB(1);
ax.vAxC = ax.vAxC-ax.vAxC(1);

[U_mm, V_mm] = meshgrid(ax.uAxC, ax.vAxC);

axP = makeimaxisfn(ImgData{2}(:,:,1,1), lambda_mm*PData(2).PDelta([1 3]));
axP.vAxB = axP.vAxB-axP.vAxB(1);
axP.vAxC = axP.vAxC-axP.vAxC(1);

[UP_mm, VP_mm] = meshgrid(axP.uAxC, axP.vAxC);

numFrames = size(ImgData{imb},4);

axI = [];
axI.uAxC_mm = linspace(axP.uAxC(1), axP.uAxC(end), size(ImgData{3},2)); 
axI.vAxC_mm = linspace(axP.vAxC(1), axP.vAxC(end), size(ImgData{3},1)); 

[UI_mm, VI_mm] = meshgrid(axI.uAxC_mm, axI.vAxC_mm);

rect = [axP.uAxC(1) axP.vAxC(1) axP.uAxC(end)-axP.uAxC(1) axP.vAxC(end)-axP.vAxC(1)]; 
imSqz = squeeze(ImgData{imb});

imCrop = [];

for f = 1:numFrames
  thisIm =  imSqz(:,:,f);
  ind = find(thisIm < eps);
  thisIm(ind) = eps;
  %imCrop(:,:,f) = imcrop(ax.uAxC, ax.vAxC, thisIm.^0.15, rect);
  imCrop(:,:,f) = interp2(U_mm, V_mm, thisIm.^0.15, UI_mm, VI_mm);
end

%UCrop_mm = imcrop(ax.uAxC, ax.vAxC, U_mm, rect);
%VCrop_mm = imcrop(ax.uAxC, ax.vAxC, V_mm, rect);
UCrop_mm = UI_mm;
VCrop_mm = VI_mm;


imToReg = [];
imToReg{1} = imCrop(:,:,1);

F = [];
FP = [];
imRegTo1 = [];
dU = []; dV = [];
for f = 1:numFrames
  imToReg{2} = imCrop(:,:,f);
  [F{f}, imRegTo1(:,:,f)] = imregdemons(imToReg{2}, ...
                                        imToReg{1}, ...
                                        'DisplayWaitbar', false);
  FP{f}(:,:,1) = interp2(UCrop_mm, VCrop_mm, F{f}(:,:,1), UP_mm, VP_mm);
  FP{f}(:,:,2) = interp2(UCrop_mm, VCrop_mm, F{f}(:,:,2), UP_mm, VP_mm);  
  dU(:,:,f) = FP{f}(:,:,1);
  dV(:,:,f) = FP{f}(:,:,2);  
end

maskExcl = (max(dV,[],3)>1);


Q = 4;

ROI_mm = [-1.6520 2.1146
           6.8133 15.6097];

% sulcul artery 
ROISmall_mm = [-0.057 0.736   
                9.495 11.17];

rectROI_mm = [ROI_mm(:,1).' ROI_mm(:,2).'-ROI_mm(:,1).'];
cImf1 = [];

figure(2)
for f = 1:numFrames
  clf
  colormap(gray)
  t = tiledlayout(2,Q,'TileSpacing','Compact','Padding','Compact');
  nexttile
  imagesc(ax.uAxB, ax.vAxB,imToReg{1}); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  nexttile
  imagesc(ax.uAxB, ax.vAxB,imCrop(:,:,f))
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b'); 
  nexttile
  imagesc(ax.uAxB, ax.vAxB,imRegTo1(:,:,f))
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  nexttile
  %  imagesc(ax.uAxB, ax.vAxB, sqrt(F{f}(:,:,1).^2+F{f}(:,:,2).^2))
    imagesc(ax.uAxB, ax.vAxB, F{f}(:,:,2))
  nexttile
  cIm1 = ImgData{2}(:,:,1,1);  
  imagesc(ax.uAxB, ax.vAxB, cIm1); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  nexttile
  cImf = ImgData{2}(:,:,1,f);  
  imagesc(ax.uAxB, ax.vAxB, cImf); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  nexttile
  cImf1(:,:,f) = imwarp(cImf, FP{f});
  imagesc(ax.uAxB, ax.vAxB, cImf1(:,:,f)); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  %  nexttile
  %imagesc(ax.uAxB, ax.vAxB, cImf1.*~maskExcl); 
  %rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
    
 
  pausede
end


makeMovie=1;
movieSpec.movieName = [vtData '/' mfile '.avi'];
vidObj = VideoWriter(movieSpec.movieName);
vidObj.FrameRate = 1;
open(vidObj);

figure(3)
for f = 2:numFrames
  clf
  t = tiledlayout(1,3,'TileSpacing','Compact','Padding','Compact');
  nexttile
  cImf = ImgData{2}(:,:,1,f);
  imagesc(ax.uAxB, ax.vAxB, cImf); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  colormap(gray)
  xlabel('mm');
  ylabel('mm');
  title('No motion correction');
  nexttile
  
  imagesc(ax.uAxB, ax.vAxB, cImf1(:,:,f)); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  colormap(gray)
  xlabel('mm');
  ylabel('mm');
  title('Motion corrected');
  
  nexttile
  imagesc(ax.uAxB, ax.vAxB, abs(cImf1(:,:,f)-cImf)); 
  rw = rectangle('Position', rectROI_mm, 'Facecolor', 'none', 'Edgecolor', 'b');
  colormap(gray)
  xlabel('mm');
  ylabel('mm');
  title('Difference');

  
   if makeMovie
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame);
  end

  pausede
end

if makeMovie
  close(vidObj);
end


return


for f = 1:numFrames
figure(1)
clf
%imagesc(log(d.ImgData{1}(:,:,1)))%imagesc(log(d.ImgData{1}(:,:,1)))
if pd == 1
  imagesc(ax.uAxB, ax.vAxB, log(ImgData{imb}(:,:,1,f)))
else
  imagesc(ax.uAxB, ax.vAxB, (ImgData{imb}(:,:,1,f)))
end

colormap gray
xlabel('mm');
ylabel('mm');
drawnow
pausede
end

save(outFile, 'ax', 'ImgData', 'lambda_mm', 'PData');
lslrt(outFile);

%figure(2)
%clf
%imshowpair(imCrop{1}(:,:,center1(3)), R2D{1}, ...
%           imCrop{2}(:,:,center2(3)), R2D{2})
           %title('Unregistered Axial slice')
           %drawnow

           %disp('Press key to start registration')
           %pausede
           %disp('Proceeding ...');


