function [rotX, rotY] = rot2d(x,y,ang)

x = x(:).';
y = y(:).';

xyMat = [ x; y];
xryrMat = [cos(ang) -sin(ang);
           sin(ang)  cos(ang)] * xyMat;
rotX = xryrMat(1,:).';
rotY = xryrMat(2,:).';


  