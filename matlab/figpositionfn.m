function figPos = figpositionfn(hImDispWin);

persistent hImDispWinSto;
  
if nargin == 0
  hImDispWin = hImDispWinSto;
else
  hImDispWinSto = hImDispWin;
end

if ~isempty(hImDispWin) & ishandle(hImDispWin)
  set(hImDispWin, 'MenuBar', 'none');
  set(hImDispWin, 'ToolBar', 'none');
end

  
figSize = figsizefn;
  
figPos = [];
screenSizeMatlab = get(0,'ScreenSize'); % this is sometimes incorrect
screenSize = getscreensizefn;
screenSize(3)=screenSize(1);
screenSize(4)=screenSize(2);  

assigned = 0;

if ~isempty(hImDispWin) & ishandle(hImDispWin)
  BOuterPosition = get(hImDispWin, 'OuterPosition');
  BPosition = get(hImDispWin, 'Position');
  setwindowbnsiconfn(hImDispWin);
else
  BPosition = [0 0 0 0];      
end

assigned = 1;
figPos.screenSize = screenSize;
figPos.screenSizeMatlab = screenSizeMatlab;

figPos.usB = [figSize.margin.S ...
              screenSizeMatlab(4)-figSize.usBSize(2)-figSize.margin.T ...
              figSize.usBSize];

offsetSD.w = figPos.usB(3) + figPos.usB(1);

figPos.usM = [offsetSD.w ...
              screenSizeMatlab(4)-figSize.usMSize(2)-figSize.margin.T ...
              figSize.usMSize];

figPos.sd1 = [offsetSD.w ...
              screenSizeMatlab(4)-figSize.sdSize(2)-figSize.margin.T ...
                figSize.sdSize];
figPos.sd2HeadGap = 29;
figPos.sd2 = [offsetSD.w ...
 screenSizeMatlab(4)-2*figSize.sdSize(2)-figSize.margin.T-figPos.sd2HeadGap ...
              figSize.sdSize];

%515.0000  409.0000
%figPos.VSX = [ 1179 431 450 620];
figPos.VSX = [ 2106 49 450 620];

figPos.xBCenter = BPosition(1)+BPosition(3)/2;

figPos.sliderLeft = figPos.xBCenter-figSize.sliderWidth/2;
figPos.sliderTop =  figPos.usB(4)-figSize.sliderHeight*1.4;
figPos.sliderTextLeft = figPos.sliderLeft;
figPos.sliderTextTop =  figPos.sliderTop-figSize.sliderTextHeight*1.2;
figPos.sliderBottom = 0;
figPos.depthSlider = [figPos.sliderLeft figPos.sliderTop ...
                    figSize.sliderWidth figSize.sliderHeight];
figPos.depthSliderTxt = [figPos.sliderTextLeft figPos.sliderTextTop ...
                           figSize.sliderTextWidth figSize.sliderTextHeight];

  %figPos.measConsole = [ 1344         282         601
  %769];
figPos.measConsole = [ 1951 597.5 615 813];
  
figPos.figSize = figSize;
figPos.waitbar = [967 554 270 56];

figPos.resultWindow = [1357.5 295 204 128];

  

  
