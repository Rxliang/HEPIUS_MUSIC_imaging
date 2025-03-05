% function ax = makedftaxis(lenAx, dbg)
% Make axis for DFT of length N points
% If dbg == 1 will check to see that 
% axis is a symmetric valid DFT axis

function [ax, foldStart] = makedftaxis(lenAx, dbg)

if nargin < 2 % If only 1 arg is given
  dbg = 0;    % set debug mode off
end

axMid = round(lenAx/2);
axStart = 0:axMid-mod(lenAx,2);
foldStart = length(axStart)+1;
if lenAx==1
  foldStart = 1;
end
axEnd = fliplr(axStart(2:(end-1+mod(lenAx,2))));
ax = [axStart axEnd];


if dbg
  F = fft(ax);
  checkreim(F);
end