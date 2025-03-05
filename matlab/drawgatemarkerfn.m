function h = ...
    drawgatemarkerfn(hnd, x1, y1, x2, y2, col, lineWid)
        

ptsx = [x1;
        x2;
        x2;
        x1];
        
ptsy = [y1;
        y1;
        y2;
        y2];

h(1) = linept2ptfn([ptsx(1); ptsy(1)], [ptsx(2); ptsy(2)], hnd); 
h(2) = linept2ptfn([ptsx(2); ptsy(2)], [ptsx(3); ptsy(3)], hnd); 
h(3) = linept2ptfn([ptsx(3); ptsy(3)], [ptsx(4); ptsy(4)], hnd); 
h(4) = linept2ptfn([ptsx(1); ptsy(1)], [ptsx(4); ptsy(4)], hnd); 

set(h, 'color', col, 'linewidth', lineWid);

