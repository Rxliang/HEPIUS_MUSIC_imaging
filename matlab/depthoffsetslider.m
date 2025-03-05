function depthoffsetslider

evParms = evalin('base', 'evParms');

persistent hl

if exist('hl', 'var') & isgraphics(hl) & isvalid(hl)
    delete(hl);
end
 
if ~evParms.state.initDepthSliderDone
    depthOffset_mm = evParms.B.depthOffset_mm;
    Resource = evalin( 'base','Resource');
    hImDispWin = Resource.DisplayWindow(1).figureHandle;
    ax = get(hImDispWin, 'CurrentAxes');
    yTickOrig = ax.YTick;
    figPos = figpositionfn;
    maxSlider_mm=50;
    slhan=uicontrol('parent', hImDispWin, ...
                    'style','slider','position', figPos.depthSlider,...
                    'min',0,'max',maxSlider_mm,'callback',@dscallbackfn,...
                    'sliderstep', [1 5]/maxSlider_mm);
    set(slhan, 'Value', 0);
    sliderStr = ['depth offset = ' num2str(depthOffset_mm) ' mm'];
    hsttext=uicontrol('parent', hImDispWin,'style','text',...
                      'position', figPos.depthSliderTxt,...
                      'visible','on', 'string', ...
                       sliderStr, 'foregroundcolor', 'k', ...
                      'fontsize', 14);
%    keyboard
    evParms.state.initDepthSliderDone=1;
    assignin('base', 'evParms', evParms);   
end

set(hsttext, 'string', sliderStr);

hl = drawlinesetaxis(yTickOrig, depthOffset_mm);

  function dscallbackfn(source,eventdata)
        
        if ~isempty(hl)
            delete(hl);
        end
        depthOffsetPre_mm = get(source,'value');
        depthOffset_mm = round(depthOffsetPre_mm);
        sliderStr = ['depth offset = ' num2str(depthOffset_mm) ' mm'];
        set(hsttext,'visible','on','string', sliderStr)
        
        if exist('hl', 'var') & isgraphics(hl) & isvalid(hl)
            delete(hl);
        end

        hl = drawlinesetaxis(yTickOrig, depthOffset_mm);
        evParms = evalin('base', 'evParms', evParms); 
        evParms.B.depthOffset_mm = depthOffset_mm;
        assignin('base', 'evParms', evParms);  

  end

  function hl = drawlinesetaxis(yTickOrig, depthOffset_mm)

    hl = line(xlim, depthOffset_mm*[1 1]);
    set(hl, 'linestyle', ':', 'linewidth',2, 'color', 'y');
        
    ytl = [];
    yt = yTickOrig;
    for i = 1:length(yt)
        ytl{i,1} = num2str(yt(i)-depthOffset_mm);
    end
    set(ax, 'YTickLabel', ytl);
  end
 
end
