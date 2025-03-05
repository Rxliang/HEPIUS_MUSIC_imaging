function nicpseqwt_roiplot(varargin)

persistent ROIHandle 
%disp('start of ROIplot');

evParms =  evalin('base','evParms');
waitHnd = findall(allchild(0), 'flat', 'Tag', 'startwaitbar');
% if multiple waitbars exist, don't meddle
if ~isempty(waitHnd) & ishandle(waitHnd) & (length(waitHnd)==1)
  waitbarExists=1;
  waitbar(0.8, waitHnd);
else
  waitbarExists=0;
end

P = evParms.P;
UI = evalin('base','UI');
UIPos = evalin('base','UIPos');
Trans = evalin('base','Trans');
PData = evalin('base','PData');
PDataIndCFI2 = evParms.ind.PDataIndCFI2;
Resource = evalin('base','Resource');
singleTpcProfile = evParms.flag.singleTpcProfile;

% change the txt for voltage slider
hv1txt = findobj('tag','hv1txt'); 
hv2txt = findobj('tag','hv2txt');
set(hv1txt,'String','B/SD Voltage');
if ~singleTpcProfile
  set(hv2txt,'String','CFI Voltage');
end

% initialize CFI Doppler state control buttons
switch evParms.state.DopState
 case 'freq'
  v1=1; v2=0;
 case 'power'
  v1=0; v2=1;
end

[exFlag, val] = existbasefn('Mcr_GuiHide');
if exFlag
  Mcr_GuiHide = val;
else
  Mcr_GuiHide = 0;
end

if ~Mcr_GuiHide
    set(findobj('tag', [evParms.UI.button.dopplermode.pos 'RadioButton1']), ...
        'value', v1);
    set(findobj('tag', [evParms.UI.button.dopplermode.pos 'RadioButton2']), ...
        'value', v2);
else
    dopplerModeHnd = findobj('Title','Doppler Mode');
    if ~isempty(dopplerModeHnd) & ishandle(dopplerModeHnd)
      set(dopplerModeHnd, 'visible', 'off');
    end    
end

if ~evParms.state.initFigsDone
  hImDispWin = Resource.DisplayWindow(1).figureHandle;
  figPos = figpositionfn(hImDispWin);
  evParms.figPos = figPos;
  set(hImDispWin, 'Position', figPos.usB); 
  ylabel(hImDispWin.CurrentAxes, 'mm');
  set(hImDispWin, 'Name', evParms.UI.USImageFigureTitle);

  VSXWin =  findobj('tag','UI'); 
  if 1
    disp('DISABLED VSX POSITIONING')
  else
    set(VSXWin, 'Position', figPos.VSX);
  end
  
  scaleToWvl = 1/evParms.gate.lambda_mm;    
  unitMm=1;
  evParms.gate.SDGateMarkerHandle = plotsdgatefn(evParms, unitMm, scaleToWvl);

  disp(['x = ' num2str(evParms.gate.x/scaleToWvl)]);
  disp(['z = ' num2str(evParms.gate.z/scaleToWvl)]);    

  evParms.state.initFigsDone = 1;   
  assignin('base', 'evParms', evParms);
        
  % put a depth offset slider on the main B mode window
  if ~isfield(evParms.flag, 'depthOffsetSliderFlag') | evParms.flag.depthOffsetSliderFlag
    depthoffsetslider; % writes modified evParms to base
  end
  
  
  % draw SD overlay
  drawsdoverlayroifn; % writes modified evParms to base
  if waitbarExists & ishandle(waitHnd) & (length(waitHnd)==1)
    waitbar(1, waitHnd);
    pause(0.25);
    close(waitHnd);
  end
end

if 0 & ~evParms.flag.doCFINow
  disp('end of ROIplot for ~doCFINow');
  return
end

if strcmp(Trans.name, 'GEC1-6D') | strcmp(Trans.name, 'GEC1-5D') 
  sector=1;
else
  sector=0;
end

if strcmp(Trans.name, 'GEM5ScD_64') | strcmp(Trans.name, 'GEM5ScD') 
  sector=0;
end

if ~sector
  PP = P(evParms.ind.PIndCFI);
  % Outline the PData region for doppler processing
  X(1) = PData(PDataIndCFI2).Region.Shape.Position(1)-PData(PDataIndCFI2).Region.Shape.width/2;
  Y(1) = PData(PDataIndCFI2).Region.Shape.Position(3);
  if Y(1) < PP.startDepth, Y(1) = PP.startDepth; end
  X(2) = X(1) + PData(PDataIndCFI2).Region.Shape.height*tan(PP.dopAngle);
  X(3) = X(2) + PData(PDataIndCFI2).Region.Shape.width;
  Y(2) = Y(1) + PData(PDataIndCFI2).Region.Shape.height;
  Y(3) = Y(2);
  X(4) = X(1) + PData(PDataIndCFI2).Region.Shape.width;
  Y(4) = Y(1);
  if PP.dopAngle > 0    
    if X(3) >  PData(PDataIndCFI2).Origin(1)+PData(PDataIndCFI2).Size(2)*PData(PDataIndCFI2).PDelta(1)
      X(5) = X(4); Y(5) = Y(4);
      X(3) = PData(PDataIndCFI2).Origin(1)+PData(PDataIndCFI2).Size(2)*PData(PDataIndCFI2).PDelta(1);
      X(4) = X(3);
      Y(4) = (256-PP.dopDispEle)/2*Trans.spacing/tan(PP.dopAngle);
    end
  else
    if X(2) < PData(PDataIndCFI2).Origin(1)
      X(5) = X(4); Y(5) = Y(4);
      X(4) = X(3); Y(4) = Y(3);
      X(3) = PData(PDataIndCFI2).Origin(1);
      X(2) = X(3);
      Y(2) = (256-PP.dopDispEle)/2*Trans.spacing/tan(-PP.dopAngle);                
    end
  end
  
  ROIX = [X X(1)];
  ROIY = [Y Y(1)];
else % sector
  % Outline the PData region for doppler processing
  if evParms.flag.usem5 |  evParms.flag.use6s 
    r1 = PData(PDataIndCFI2).Region.Shape.z;
    r2 = PData(PDataIndCFI2).Region.Shape.r;
    % approximation of sectorFT
  else  
    r1 = PData(PDataIndCFI2).Region.Shape.r1;
    r2 = PData(PDataIndCFI2).Region.Shape.r2;
  end
      
  angle = PData(PDataIndCFI2).Region.Shape.angle; % in radian
  steerAngle = P(evParms.rcvOut.PIndCFI).dopAngle; % in radian
  theta = pi/2-angle/2-steerAngle; % start angle
  circlePts = 180; % points for whole circle
  n = abs(ceil(circlePts*angle/(2*pi)));
  curve_theta = theta + (angle*(0:n)'/n);
  curve1_R = r1*ones(n+1,1);
  curve2_R = r2*ones(n+1,1);

  [curve1_X,curve1_Y] = pol2cart(curve_theta,curve1_R);
  [curve2_X,curve2_Y] = pol2cart(curve_theta,curve2_R);

  % correct ROI points if curve2 is outside of the boundary
  xmin = Resource.DisplayWindow(1).ReferencePt(1);
  xmax = xmin + Resource.DisplayWindow(1).Position(3)*Resource.DisplayWindow(1).pdelta;
  if min(curve2_X) < xmin
    X1 = curve1_X(end); Y1 = curve1_Y(end);
    X2 = curve2_X(end); Y2 = curve2_Y(end);
    X = xmin; Y = Y1 + (X-X1)*(Y2-Y1)/(X2-X1);
    curve2_X(curve2_X<xmin) = xmin;
    curve2_X(end+1) = X;
    curve2_Y(end+1) = Y;
  elseif max(curve2_X) > xmax
    X1 = curve1_X(1); Y1 = curve1_Y(1);
    X2 = curve2_X(1); Y2 = curve2_Y(1);
    X = xmax; Y = Y1 + (X-X1)*(Y2-Y1)/(X2-X1);
    curve2_X(curve2_X>xmax) = xmax;
    curve2_X = [X; curve2_X];
    curve2_Y = [Y; curve2_Y];
  end
  keyboard
  ROIX = [curve1_X; flipud(curve2_X); curve1_X(1)]+PData(PDataIndCFI2).Region.Shape.Position(1); % shift X
  ROIY = [curve1_Y; flipud(curve2_Y); curve1_Y(1)]+PData(PDataIndCFI2).Region.Shape.Position(3); % shift Y
end %

% scale to mm if the axes unit of the displaywindow is mm
scaleToMm = (Resource.Parameters.speedOfSound/1000)/Trans.frequency;
if isfield(Resource.DisplayWindow(1),'AxesUnits')&&~isempty(Resource.DisplayWindow(1).AxesUnits)
    if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
        ROIX = ROIX * scaleToMm; ROIY = ROIY * scaleToMm;
    end
end
bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');
if ishandle(bmodeFigHandle)
    if  isempty(ROIHandle) || ~ishandle(ROIHandle)
        figure(bmodeFigHandle), hold on,
        ROIHandle = plot(ROIX,ROIY,'Color',[0.7 0.7 0.7],...
			 'LineWidth',1.5,'LineStyle','--','Visible','off');
        if get(UI(6).handle(2),'Value') == 1 % 'on' is selected
            set(ROIHandle,'Visible','on');
        end
        assignin('base','ROIHandle',ROIHandle);
    else
        set(ROIHandle,'XData',ROIX);
        set(ROIHandle,'YData',ROIY);
    end
end

%disp('end of ROIplot');

%return

end

