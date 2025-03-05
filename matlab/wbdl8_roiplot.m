function wbdl8_roiplot(varargin)

persistent drawROI

% change the txt for voltage slider
hv1txt = findobj('tag','hv1txt');
hv2txt = findobj('tag','hv2txt');
set(hv1txt,'String','Bmode Voltage');
set(hv2txt,'String','Doppler Voltage');

% drawRegionOutline(WinNum,PDataNum,RegionNum) creates an outline from
% PData(PDataNum).Region(RegionNum) on displayWindow(WinNum) with default color - white.
if isempty(drawROI)
    drawROI = 1;
    evalin('base','hROI = drawRegionOutline(1,3,1);')
else
    evalin('base','drawRegionOutline(hROI,1,3,1);')
end

end
