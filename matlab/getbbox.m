function bbox = getbbox(x,z, wbAxes, xzDimROI)
  if nargin < 4
    xzDimROI = [1 1];
  end
  
  axun = get(wbAxes,'Units'); % save axis units
  children = get(wbAxes, 'children');
  types = get(children, 'Type');
  ind = find(strcmp(types, 'image'));
  child = children(ind);
  imSize = size(child.CData);
  set(wbAxes,'Units','normalized'); % make axes units normalized
  axpos = get(wbAxes,'Position'); % get axes position
  axlim = axis(wbAxes); % get the axis limits [xlim ylim (zlim)]
  axwidth = diff(axlim(1:2));
  axheight = diff(axlim(3:4));
  abox = [x,z,(xzDimROI(1)*axwidth)/imSize(2),(xzDimROI(2)*axheight)/imSize(1)];
  % Transform from data space coordinates to normalized figure coordinates
  bbox(1) = (abox(1)-axlim(1))/axwidth * axpos(3) + axpos(1);
  bbox(2) = (axlim(4)-abox(2))/axheight * axpos(4) + axpos(2);
  bbox(3) = abox(3) * axpos(3)/axwidth;
  bbox(4) = abox(4) * axpos(4)/axheight;
  set(wbAxes,'Units',axun); % restore axis units
end