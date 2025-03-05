function setwindowbnsiconfn(hnd)

if ~isempty(hnd) & ishandle(hnd)
  warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
  vtHome = getenv('BNS_HOME');
  iconFile = fullfile(vtHome, 'matlab', 'bnsicon_crop.gif');
  jFrame=get(hnd, 'javaframe');
  jicon=javax.swing.ImageIcon(iconFile);
  jFrame.setFigureIcon(jicon);
end
