function dispExtern(img)
persistent lasttic figH Trans PData axesHandle lastsize imgH

if isempty(lastsize),PData = evalin('base', 'PData'); lastsize = PData(1).Size(1); end

% Initialize
if isempty(lasttic)
    P = evalin('base','P');
    Trans = evalin('base','Trans');
    figH = figure('Position',[230, 240, 900, 900*(P.endDepth - P.startDepth)/(128*Trans.spacing)],'Color','k');
    set(figH, 'MenuBar', 'none');
    set(figH, 'ToolBar', 'none');
    pos = [0.1 .1 0.86 0.86]; % left bottom width height
    axesHandle = subplot('Position', pos);
    imagesc(flipud(20*log10(abs(squeeze(img)))) - 77,'parent',axesHandle,[0 56]);
    imgH = get(axesHandle, 'Children');
    set(figH, 'InvertHardcopy', 'off');
    axis equal; axis tight; colormap gray;
    set(figH,'name','Custom Live Display')
end

if (isempty(lasttic) || PData(1).Size(1) ~= lastsize)
    PData = evalin('base', 'PData');
    P = evalin('base','P');
    title(
    lasttic = tic;
    daspect(axesHandle, [PData(1).PDelta(3) PData(1).PDelta(1) 1]);
    set(axesHandle,'NextPlot','replaceChildren')
    axis(axesHandle,[1 round((PData(1).Size(2))) 1 (PData(1).Size(1))]) % width height
    set(axesHandle,'XTick',linspace(1,(PData(1).Size(2))- 1,9))
    set(axesHandle,'XTickLabel', num2cell(round(linspace(-PData(1).Size(2)/2*PData(1).PDelta(1),PData(1).Size(2)/2*PData(1).PDelta(1),9)*Trans.spacingMm,1)))
    set(axesHandle,'YTick',linspace(1,(PData(1).Size(1)) - 1,10))
    set(axesHandle,'YTickLabel', num2cell(round(linspace(P.startDepth,P.endDepth,10)*1540/1000/(Trans.frequency),1)))
end

set(imgH,'CData',20*log10(abs(squeeze(img)))- 77);

end

%-- Initialization --%
%if isempty(lasttic)
%    P = evalin('base','P');
%    Trans = evalin('base','Trans');
%    figH = figure('Position',[230, 240, 900, 900*(P.endDepth - P.startDepth)/(128*Trans.spacing)],'Color','k');
%    set(figH, 'MenuBar', 'none');
%    set(figH, 'ToolBar', 'none');
%    pos = [0.1 .1 0.86 0.86]; % left bottom width height
%    axesHandle = subplot('Position', pos);
%    imagesc(flipud(20*log10(abs(squeeze(complex(IData, QData)))))- 115,'parent',axesHandle,[0 50]);
%    imgH = get(axesHandle, 'Children');
%    set(figH, 'InvertHardcopy', 'off');
%    axis equal; axis tight; colormap gray;
%    set(figH,'name','Custom Live Display')
%    PData = evalin('base', 'PData');
%    imagedim=PData(1).Size;
%    Z=PData(1).PDelta(3)*(0:imagedim(1)-1);
%    X=PData(1).PDelta(1)*(0:imagedim(2)-1);
%    dopPRF = evalin('base','dopPRF');
%    TX = evalin('base','TX');
%    nPM = evalin('base','nPM');
%    nAngles = evalin('base', 'nAngles');
%    myangles = [];
%    for x = 1:nPM:nPM*nAngles
%        myangles = [myangles TX(x).Steer(1)*180/pi];
%    end
%    sequence = evalin('base', 'sequence');
%    t1 = title([sequence ', Acquisition PRF: ' num2str(dopPRF/1000) ' kHz' ...
%        '       Angles: ' num2str(round(myangles)) ' degrees']);
%    t1.Color = [1 1 1];
%    set(gca,'XTick',linspace(0,X(end),7))
%    set(gca,'XTickLabel', num2cell(round(linspace(-PData(2).Size(2)/2*PData(2).PDelta(1),PData(2).Size(2)/2*PData(2).PDelta(1),7)*Trans.spacingMm,1)))
%    set(gca,'YTick',linspace(0,Z(end),5))
%    set(gca,'YTickLabel', num2cell(round(linspace(P.startDepth,P.endDepth,5)*1540/1000/(Trans.frequency),1)))
%    set(gca, 'XColor','w','YColor','w');
%    xlabel('mm','Color','w'); ylabel('mm','Color','w');
%    daspect(axesHandle, [PData(1).PDelta(3) PData(1).PDelta(1) 1]);
%    caxis([0 50]);
%    h = colorbar('Color','w'); ylabel(h, 'dB', 'Color', 'w');
%    assignin('base','figH',figH);
%end
%
%
% %-- if updated --%
%if (isempty(lasttic) || PData(1).Size(1) ~= lastsize)
%    PData = evalin('base', 'PData');
%    P = evalin('base','P');
%    lasttic = tic;
%    daspect(axesHandle, [PData(1).PDelta(3) PData(1).PDelta(1) 1]);
%    set(axesHandle,'NextPlot','replaceChildren')
%    axis(axesHandle,[1 round((PData(1).Size(2))) 1 (PData(1).Size(1))]) % width height
%    set(axesHandle,'XTick',linspace(1,(PData(1).Size(2))- 1,9))
%    set(axesHandle,'XTickLabel', num2cell(round(linspace(-PData(1).Size(2)/2*PData(1).PDelta(1),PData(1).Size(2)/2*PData(1).PDelta(1),9)*Trans.spacingMm,1)))
%    set(axesHandle,'YTick',linspace(1,(PData(1).Size(1)) - 1,10))
%    set(axesHandle,'YTickLabel', num2cell(round(linspace(P.startDepth,P.endDepth,10)*1540/1000/(Trans.frequency),1)))
%end
%
% %---Every time-----------------------------------------------%
% % Jonah's faster display
%set(imgH,'CData',20*log10(abs(squeeze(complex(IData, QData))))- 115);
