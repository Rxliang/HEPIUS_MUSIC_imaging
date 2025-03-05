function drawsdoverlayroifn

evParms = evalin('base', 'evParms');

X_mm(1) = evParms.gate.SDLarge.xSDOffset_mm- ...
          evParms.gate.SDLarge.xSDWidth_mm/2;

Y_mm(1) = evParms.gate.SDLarge.zSDStart_mm;

X_mm(2) = evParms.gate.SDLarge.xSDOffset_mm + ...
          evParms.gate.SDLarge.xSDWidth_mm/2;

Y_mm(2) = evParms.gate.SDLarge.zSDStart_mm + ...
          evParms.gate.SDLarge.zSDHeight_mm;

ROIX = [X_mm(1) X_mm(1) X_mm(2) X_mm(2) X_mm(1)];
ROIY = [Y_mm(1) Y_mm(2) Y_mm(2) Y_mm(1) Y_mm(1)];   

bmodeFigHandle = evalin('base','Resource.DisplayWindow(1).figureHandle');

if ~isfield(evParms.UI, 'SDROIHandle')
  evParms.UI.SDROIHandle=[];
end

if ishandle(bmodeFigHandle)
    if  isempty(evParms.UI.SDROIHandle) || ~ ...
          ishandle(evParms.UI.SDROIHandle)
        figure(bmodeFigHandle), hold on,
        evParms.UI.SDROIHandle = plot(ROIX,ROIY,'Color',[0.7 0.7 0.7],...
                           'LineWidth',1.5,'LineStyle','--','Visible','off');
        assignin('base', 'evParms', evParms);
    else
        set(evParms.UI.SDROIHandle,'XData',ROIX);
        set(evParms.UI.SDROIHandle,'YData',ROIY);
    end
    if evParms.largeSDParms.SDLargeDisplayOverlay
      set(evParms.UI.SDROIHandle,'Visible','on');
    else
      set(evParms.UI.SDROIHandle,'Visible','off');
    end
end


