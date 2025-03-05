function evParms = setoverlaysdpriorityfn(evParms)

%fix overlay size and UI first
hnd=findguihandlefn(evParms.UI.UCCTitle);
                
if ~isempty(hnd) & ishandle(hnd)
    txt = ...
        hnd.RunningAppInstance.ScanpriorityButtonGroup.SelectedObject.Text;
    if ~evParms.flag.doWT
      if strcmp(txt, 'Imaging') | strcmp(txt, 'Wall track')
        hnd.RunningAppInstance.ScanpriorityButtonGroup.SelectedObject.Value=0;
      end 
    else
      if strcmp(txt, 'Imaging') | strcmp(txt, 'Spectrogram')
        hnd.RunningAppInstance.ScanpriorityButtonGroup.SelectedObject.Value=2;
      end 
    end
          
    hnd.RunningAppInstance.BmodeimagePanel.Visible = 'Off';
    hnd.RunningAppInstance.RestoredefaultsBmodeButton.Enable = 'Off';    
else
  return
end
          
[gateOut, changeFlag] = autoadjustoverlayfn(evParms.gate);

if ~isempty(changeFlag) 
    evParms.gate = gateOut;
    assignin('base', 'evParms', evParms);
else
    changeFlag = [0 0 0 0];
end

dum=1;
if ~isempty(find(changeFlag))            
    %gateOut.SDLarge.xSDOffset_mm
    %gateOut.SDLarge.xSDWidth_mm
    %gateOut.SDLarge.zSDStart_mm
    %gateOut.SDLarge.zSDHeight_mm
            
  hnd.RunningAppInstance.LateraloffsetmmSlider.Value ...
      = gateOut.SDLarge.xSDOffset_mm; ...              
  hnd.RunningAppInstance.WidthmmSlider.Value ...
      = gateOut.SDLarge.xSDWidth_mm; ...
  hnd.RunningAppInstance.StartdepthmmSlider.Value=...
      gateOut.SDLarge.zSDStart_mm;                      
  hnd.RunningAppInstance.HeightmmSlider.Value=...
      gateOut.SDLarge.zSDHeight_mm;  
                
  changeParms.changeLateralOffset_mm = ...
      gateOut.SDLarge.xSDOffset_mm;
  changeParms.changeLargeSDWidth_mm = ...
      gateOut.SDLarge.xSDWidth_mm;               
  changeParms.largeSDHeight_mm = ...
      gateOut.SDLarge.zSDHeight_mm;
  changeParms.largeSDZStart_mm = gateOut.SDLarge.zSDStart_mm;               
  changeParms.flag.changeLargeSDWidth=1;
  changeParms.flag.changeLargeSDLateralOffset=1;
  changeParms.flag.changeLargeSDHeight=1;
  changeParms.flag.changeLargeSDZStart=1;
  evParms = updatepdatactrlfn(changeParms); % assigned
                                            % to base already
  drawsdoverlayroifn; 
end % change of x and z parameters                  
