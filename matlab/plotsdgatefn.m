function markerHndl = plotsdgatefn(evParms, unitMm, scaleToWvl)

   %markerColor = [1   0.2 .2;
%               0.2   1 .2];

% row 1: blue for proximal gate    
% row 2: green for distal gate    

switch evParms.state.CDIState 
  case {1,3}
    markerColor = evParms.SD.markerColorFlow;    
  case {2,4}
    markerColor = evParms.SD.markerColorPower;
end


wbFig = evalin('base','Resource.DisplayWindow(1).figureHandle');
wbAxes = get(wbFig,'CurrentAxes');

for q = 1:evParms.gate.numSD
   x_mm = evParms.gate.x(q)*evParms.gate.lambda_mm;
   z_mm = evParms.gate.z(q)*evParms.gate.lambda_mm;
   
   rx1 = x_mm-evParms.gate.xSVWidth_mm/2;
   rx2 = x_mm+evParms.gate.xSVWidth_mm/2;

   rz1 = z_mm-evParms.gate.zSVWidth_mm/2;
   rz2 = z_mm+evParms.gate.zSVWidth_mm/2;   
   col =  markerColor(q,:);
   lineWid = 2;
   markerHndl{q} = ...
    drawgatemarkerfn(wbAxes, rx1,rz1,rx2,rz2, col, lineWid);

end
