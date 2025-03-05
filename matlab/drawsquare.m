function [h,ht,rx,ry] = ...
    drawsquare(x1,y1,x2,y2, label, col, fillIn, angDeg)

h=[];
ht=[];

if nargin < 6
  col = 'k';
end
if nargin < 7
  fillIn = 0;
end

if nargin < 8
  angDeg=0;
end

ptsx = [x1;
        x2;
        x2;
        x1];
        
ptsy = [y1;
        y1;
        y2;
        y2];
        
[rx, ry] = rot2d(ptsx,ptsy, angDeg/180*pi);
      
if strcmp(col, 'none')
  return
end

h(1) = linept2pt([rx(1); ry(1)], [rx(2); ry(2)]); 
h(2) = linept2pt([rx(2); ry(2)], [rx(3); ry(3)]); 
h(3) = linept2pt([rx(3); ry(3)], [rx(4); ry(4)]); 
h(4) = linept2pt([rx(1); ry(1)], [rx(4); ry(4)]); 

set(h, 'color', col);

%h(1) = line([x1 x2], [y1 y1],'color',col);
%h(2) = line([x2 x2], [y1 y2],'color',col);
%h(3) = line([x2 x1], [y2 y2],'color',col);
%h(4) = line([x1 x1], [y2 y1],'color',col);

if fillIn
  X = [x1 x2 x2 x1];
  Y = [y1 y1 y2 y2];
  hold on;
  fill(X,Y,col);
  hold off;
end  

if nargin > 4 & ~isempty(label)
  midPtX = ((x1+x2)/2)- .25*abs(x2-x1)*length(label)/3;
  midPtY = (y1+y2)/2-.03*abs(x2-x1);
  ht = text(midPtX, midPtY, label, 'fontsize', abs(x2-x1), 'color',col)
end
