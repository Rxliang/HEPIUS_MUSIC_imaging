clear files spectralDoppler_4p9p2_saveiq
global elypHandlGlobal 
interactiveMode = 0;

% run offlinesdlarge first

hFigMain = figure(figNoIm)

clf
tld = tiledlayout(2,1,'TileSpacing','Compact','Padding','Compact');
hTile = [];
hTile.tile = tld;
hTile.fig = hFigMain;

colormap(gray)
set(gcf, 'WindowButtonDownFcn', @wbdCallback);
set(gcf, 'position', [35  560  852 497]);
nexttile

%set(gcf, 'position', [10.0000  850.8571  852 497]);

hTile.imBCDI = imagesc(iSet{1}.xAx_mm, iSet{1}.zAx_mm, imCompNorm(:,:,:,1));
xlim(xlm); ylim(zlm);

%imagesc(iqSet.xAx_mm, iqSet.zAx_mm, imBCDI);
xlabel('x (mm)');
ylabel('z (mm)');
axis equal

ind = find(flowImNorm == maxall(flowImNorm));
[r,c] = ind2ij(ind, size(flowImNorm,1));
r=r(1);
c=c(1);
%r=r-20;
x = iqSet.xAx_mm(c); % axes of IQ inset
z = iqSet.zAx_mm(r);
wbAxes = get(gcf,'CurrentAxes');
bbox = getbbox(x,z, wbAxes);
%elypHndlGlobal = annotation('rectangle',bbox,'EdgeColor','y','LineWidth',2.0);

[dum, rLim] = findclosestinvec(iqSet.zAx_mm, xzRect{s}(:,2));
[dum, cLim] = findclosestinvec(iqSet.xAx_mm, xzRect{s}(:,1));

nPulses = size(iqSet.IQMat,3);

if length(rLim) > 1
  rRange = rLim(1):rLim(2);
else
  rRange = rLim;
end

if length(cLim) > 1
  cRange = cLim(1):cLim(2);
else
  cRange = cLim;
end

indFlowImNorm = reshape(1:prod(size(flowImNorm)), size(flowImNorm));

cDopROI = flowImNorm(rRange, cRange);
cDopROIInd = indFlowImNorm(rRange, cRange);

[cDopSort, sortInd] = sort(cDopROI(:));

% how many pixels to keep in measured IQ region
keepNTop_used = min(length(sortInd), keepNTop);
topInd = cDopROIInd(sortInd(end-keepNTop_used+1:end));

setInd = s;
save(outMatFile, 'topInd', 'xzRect', 'sgSel', 'inPath', 'music_sets', 'outMatFolderFile', 'setInd', 'large_SD_dim');
lslrt(outMatFile)

if 0
IQTemp = iqSet.IQMat;
len2 = prod(size(IQTemp,[1 2]));
len3 = prod(size(IQTemp));
descale = 2;
for q = 1:length(topInd)
  ind = topInd(q):len2:len3;
  thisSig = IQTemp(ind);
  iqdr = wdenoise(real(thisSig), descale);
  iqdi = wdenoise(imag(thisSig), descale);        
  %iqSet.IQMat(ind) = iqdr+j*iqdi;
  figure(10)
  ax = 1:length(thisSig);
  subplot(1,3,1)
  plot(ax, real(thisSig), ax, imag(thisSig));
  subplot(1,3,2)
  plot(ax, iqdr, ax, iqdi)
  subplot(1,3,3)
  plot(ax, iqdr-real(thisSig), ax, iqdi-imag(thisSig))  
  %pausede
end

end

iqSet.evParms.SD.SDDynRangeDB = 3;
iqSet.evParms.SD.despeckle = 0.4;
iqSet.evParms.SD.noisePersist = 0.95/2;

iqSet.evParms.SD.narrowbandInd = topInd; % use these

if ~interactiveMode
  iqSet.evParms.SD.narrowbandFlag = 1; % use narrowbanding
  iqSet.evParms.SD.narrowbandIndFlag = 1; % using provided indices
else
  iqSet.evParms.SD.narrowbandFlag = 0; % don't use narrowbanding  
end

% if these green dots don't show up at true maximum flow positions,
% it is likely xzRect is not set appropriately for this dataset
% and search for max flow power is taking place in wrong ROI

hold on
iqSet.evParms.SD.hplotTopInd = plot(iqSet.X_mm(topInd), iqSet.Z_mm(topInd), 'g.');
%rectangle('Position', [iqSet.Z_mm(topInd), iqSet.X_mm(topInd) 1 1]);
hold off
 
sdOut = [];
sdOut.Sxx = [];
sdOut.IQSlideBlocks = [];
sdOut.SDop = [];
sdOut.sgp = [];
sdOut.saveFlag = 1;
sdOut.offset=1;

init = 1;
while 1
    %[r(1) c]
    %r = round(r);
    %c = round(c);
    %  r = rRange;
    %c = cRange;
  if ~interactiveMode
    % feed entire IQMat to SpectralDoppler, where narrowbanding may or may not be performed
    r = 1:size(iqSet.IQMat,1);
    c = 1:size(iqSet.IQMat,2);
  end  

  disp([mfile ': interactiveMode = ' num2str(interactiveMode)]);

  disp(['Left button: Select individual pixel. Right button: Show mean signal for green pixels']);
  
  % retu n
  stopFlag = 0;
  restartFlag = 0;
  %  [r,c]
  sdOut.sgp = offlinesd4p9p2fn(iqSet, iqSet.IQMat, r, c, init);
  init=0;
 % break
end

function wbdCallback(varargin)
% Window ButtonDown callback for sample volume positioning.

    persistent init wbFig wbAxes elypHndl 

    if exist('elypHandlGlobal', 'var')
        delete(elypHandlGlobal);
    end
    
    if 0 & isempty(init)
        init = 1;
        
        elypHndlExists = evalin('base', 'exist(''elypHndl'', ''var'')');
        if elypHndlExists
          elypHndl = evalin('base', 'delete(elypHndl)');
          %          delete(elypHndl);
        else    
          elypHndl = 0;
        end
    
    end

    % Get pos of mouse click and ensure the unit is wavelength for spectral doppler
    hObject = varargin{1};
    wbFig = evalin('base', 'figNoIm');
    hh = get(wbFig);
    CDataUpper = hh.Children(1).Children(1).Children.CData;
    wbAxes = get(wbFig,'CurrentAxes');
    cp = get(wbAxes,'CurrentPoint');
    x = cp(1,1);
    z = cp(1,2);
    
    % Change Event Sequence if left mouse button is pressed
    if strcmp(hObject.SelectionType,'normal')
       % inspect single pixels
       evalin('base', 'evParms.interactiveMode = 1;');
       evalin('base', 'evParms.SD.narrowbandFlag = 0;');        
       evalin('base', 'iqSet.evParms.SD.narrowbandFlag = 0;');
        
        %        assignin('base', 'stopFlag', 1);
       xAx_mm = evalin('base', 'iqSet.xAx_mm');
       zAx_mm = evalin('base', 'iqSet.zAx_mm');       
       [~, r] = findclosestinvec(zAx_mm, z);
       [~, c] = findclosestinvec(xAx_mm, x);
       assignin('base', 'r', r(1));
       assignin('base', 'c', c(1));
       % this will make caller send single point IQmat to this function
       assignin('base', 'interactiveMode', 1);
       
       if elypHndl ~= 0
          delete(elypHndl); % delete previous annotation circle
       end
       x=xAx_mm(c);
       z=zAx_mm(r);
       disp(['r = ' num2str(r) ', c = ' num2str(c) ', x = ' num2str(x) ', z = ' num2str(z)])
       bbox = getbbox(x,z, wbAxes);

       elypHndl = annotation('rectangle',bbox,'EdgeColor','y','LineWidth',2.0);
        
    elseif strcmp(hObject.SelectionType,'alt')
        % delete annotation circle
        if elypHndl ~= 0
            delete(elypHndl);
            elypHndl = 0;
        end
        assignin('base', 'stopFlag', 1);
        % go back to non-interactive mode
        assignin('base', 'interactiveMode', 0);
    elseif strcmp(hObject.SelectionType,'extend')
      disp('middle');
      assignin('base', 'restartFlag', 1);
      
    end
end
