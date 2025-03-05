function nicpseq_wdbCallback(hObject,eventdata)
persistent init wbFig wbAxes 

evParms = evalin('base', 'evParms');
Process = evalin('base', 'Process');

if isempty(init)
    wbFig = evalin('base','Resource.DisplayWindow(1).figureHandle');
    wbAxes = get(wbFig,'CurrentAxes');
    init = 1;   
end

% if left mouse button
%get(hObject,'SelectionType')
clickEvent = get(hObject,'SelectionType');
eventStrSV1 = 'normal';
eventStrSV2 = 'alt';
eventStrToggleMode = 'extend';

selSV1 = strcmp(clickEvent, eventStrSV1);
selSV2 = strcmp(clickEvent, eventStrSV2) & (evParms.gate.numSD > 1);

if selSV1 | selSV2
   for q = 1:evParms.gate.numSD
     ind = find(ishandle(evParms.gate.SDGateMarkerHandle{q}));
     if ~isempty(ind)
       delete(evParms.gate.SDGateMarkerHandle{q}(ind));
     else
       disp(['evParms.gate.SDGateMarkerHandle{' ...
             num2str(q) '} empty']);
     end
   end
   
   cp = get(wbAxes,'CurrentPoint');
   % these are the center of the gate
   x_mm = cp(1,1);
   z_mm = cp(1,2);

   disp(['selSV1: ' num2str(selSV1) ', selSV2: ' num2str(selSV2)]);
   disp(['x_mm: ' num2str(x_mm) ', z_mm: ' num2str(z_mm)]);   
   
   Trans = evalin('base','Trans');
   Resource = evalin('base','Resource');
  
   if isfield(Resource.DisplayWindow(1),'AxesUnits')&&...
      ~isempty(Resource.DisplayWindow(1).AxesUnits)
       if strcmp(Resource.DisplayWindow(1).AxesUnits,'mm');
           scaleToWvl = evParms.Trans.frequency/...
               (Resource.Parameters.speedOfSound/1000);   
           wl2mm = (Resource.Parameters.speedOfSound/1000)/evParms.Trans.frequency;
           cp = cp*scaleToWvl;
           unitMm = 1;
       end
   else
       unitMm = 0;
   end   
   
   x = cp(1,1); % in wavelengths
   z = cp(1,2); % in wavelengths
   %[x z]   
  
   % Check for sample volume outside of measurement region.
   %ax = evalin('base', 'ax');
   
   PData = evalin('base', 'PData');   
   
   % bring in old values
   xSV_wvl = evParms.gate.x;
   zSV_wvl = evParms.gate.z;

   if selSV1
     selInd=1;
     PDataIndSDThis = evParms.gate.PDataIndSDStart;
   else
     selInd=2;
     PDataIndSDThis = evParms.gate.PDataIndSDStart+1;
   end
 
   ax = gca(Resource.DisplayWindow(1).figureHandle);
   xLim_wvl = ax.XLim*scaleToWvl;
   yLim_wvl = ax.YLim*scaleToWvl;   
  
   % check click is in bounds of figure
   if (x <= xLim_wvl(2)) & (x >= xLim_wvl(1)) & ...
      (z <= yLim_wvl(2)) & (z >= yLim_wvl(1)) 
     xSV_wvl(selInd) = x;
     zSV_wvl(selInd) = z;
   end

   % adujsting PData for these Doppler gates
   [PData, evParms] = setpdatasdfn(PData, xSV_wvl, zSV_wvl, evParms);
  
   gateSave = evParms.gate;
   gateSave.unitMm = unitMm;
   gateSave. scaleToWvl = scaleToWvl;
   save(evParms.gate.gatesFile, 'gateSave'); 
  
   evParms.gate.SDGateMarkerHandle = plotsdgatefn(evParms, unitMm, scaleToWvl);
   
   assignin('base','evParms', evParms); % change gate
   assignin('base','PData', PData);

   assignin('base', 'PDataIndSD', PDataIndSDThis);

   evalin('base','PData(PDataIndSD).Region=computeRegions(PData(PDataIndSD));');   

   % seems this restarts seq, runs ROIPlot again
   
   Control = [];
   Control(1).Command = 'update&Run';
   Control(1).Parameters = {'PData'};
   assignin('base','Control', Control);
   
end % left mouse

% middle mouse button
if  strcmp(get(hObject,'SelectionType'), eventStrToggleMode)
  
  Control = [];
  Resource = evalin('base','Resource');
  hDisplay = Resource.DisplayWindow(1).figureHandle;
  currentMap = get(hDisplay,'Colormap');
    
  % toggle the current mode
  switch evParms.flag.currentModeBCFI
   case 1
    evParms.flag.currentModeBCFI=0;
    evIndStart = evParms.ev.evIndStartB;
   case 0
    evParms.flag.currentModeBCFI=1;
    evIndStart = evParms.ev.evIndStartBCFI;
  end
  assignin('base','evParms', evParms);
  c=1;  
  if 0 % this does not help grayscale shift issue
    newMap = grayscaleCPAmap;
    newMap(1:128,:) = currentMap(1:128,:);
    
    % this appears to have no effect, needed to use lower map for
    % all B
    Control(c).Command = 'set&Run';
    Control(c).Parameters = {'DisplayWindow',1,'colormap',newMap};
    
    c=c+1;
    Control(c).Command = 'set&Run';
    Control(c).Parameters = {'ImageBuffer',1,'lastFrame',0};
    c=c+1;
  end
    

  Control(c).Command = 'set&Run';
  Control(c).Parameters = {'Parameters',1,'startEvent', evIndStart};
  evalin('base',['Resource.Parameters.startEvent = ' ...
                 num2str(evIndStart) ';']);
  assignin('base','Control',Control);
    
end

% if shift-l or r-mouse button ...
if strcmp(get(hObject,'SelectionType'),'extend')
  % toggle SDOnlyInPlace mode
  
  if evParms.state.SDOnlyInPlace
    evParms.state.SDOnlyInPlace=0;
    disp('Toggling SDOnlyInPlace off');
    dopPRFNom = evParms.ev.dopPRFImagingPriority;
  else
    evParms.state.SDOnlyInPlace=1;
    disp('Toggling SDOnlyInPlace on');
    dopPRFNom = evParms.ev.dopPRFInPlace;
  end
  assignin('base', 'evParms', evParms);
  UIValue = dopPRFNom;
  %MARK

  if evParms.flag.debugSeqSwitchMode
    changeprffn_load(UIValue);
  else
    changeprffn(UIValue);
  end
  
end
return


end

