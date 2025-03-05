% https://stackoverflow.com/questions/47431214/how-to-get-the-handles-of-all-open-app-designer-apps-uifigures

function hnd = findguihandlefn(name, closeFlag, thisAppInstanceUIHandle)

if nargin < 2 
  closeFlag = 0;
end

hnd=[];
hFigs = findall(groot, 'Type', 'figure');

% Turn off warnings:
ws(2) = warning('query','MATLAB:structOnObject');
ws(1) = warning('query','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
for indW = 1:numel(ws)
  warning('off', ws(indW).identifier);
end
% Call function:
tf = isUIFigure(hFigs);
% Restore the warnings' state:
warning(ws);

indGUI = find(tf);

if isempty(indGUI)
  return;
end

if closeFlag
  for i = 1:length(indGUI)
    hnd = hFigs(indGUI(i));
    if ishandle(hnd) & hnd~=thisAppInstanceUIHandle
      close(hnd);
    end
  end
  return
end

if nargin==0
  if length(indGUI) > 1
    warning('More than one GUI!');
  end
  hnd = hFigs(indGUI);
  return
end

hFigsGUIName = {hFigs(indGUI).Name};

if length(hFigsGUIName)>1
  warning([mfilename ': multiple UI windows ' ...
          'present! Multiple apps running!']);
end

indMatch = strcmp(hFigsGUIName, name);

if ~isempty(indMatch)
  if length(indMatch) > 1 
    warning([mfilename ': multiple UI windows  "' ...
             hFigs(indGUI(1)).Name '" present!'])
    indMatch = indMatch(end);
  end

  hnd = hFigs(indGUI(indMatch));
end

return

function tf = isUIFigure(hFigList)
  tf = arrayfun(@(x)isstruct(struct(x).ControllerInfo), hFigList);
