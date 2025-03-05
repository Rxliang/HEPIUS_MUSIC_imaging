function vecStr = vec2commadelim(vec, crInterval, marker, fmt)

if nargin < 2  | isempty(crInterval)
  crInterval = 0;
end

if nargin < 3 | isempty(marker)
  marker = ', '; % for compatibility
end

if nargin < 4
  customFormat = 0;
else
  customFormat = 1;
end

vecStr = [];

for i = 1:length(vec)
  if customFormat
    vecStr = [vecStr num2str(vec(i), fmt) marker];
  else
    vecStr = [vecStr num2str(vec(i)) marker];
  end
  
  sp=0;
  if ~mod(i,crInterval)
    vecStr = [vecStr '\\'];
    sp = 2;
  end
end

vecStr = vecStr(1:end-length(marker)-sp);
