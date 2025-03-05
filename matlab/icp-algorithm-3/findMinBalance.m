function [ICP, Factor, polynom_data, Q] = findMinBalance(Pressures0, Pressures, input, pOrder)

% Function from TCD3D to find the pressure corresponding to the minimum of
% input
% settings
threshold = 10; %0.045; 
%threshold = 0.045; 
pBias = 0.5;
% saturate
if threshold > min(input)
    inx = find(input > threshold);
else
    inx = find(input > mean(input));
    threshold = mean(input);
end

if isempty(inx)==0
    input(inx) = ones(1,length(inx))*threshold;
    pOrder = 3;
end
if input(2)==threshold
    input(1) = threshold;
end

ICP = 0;
Factor = 0;

if length(Pressures)<=2
    pOrder = 1;
end

[PP,SS] = polyfit(Pressures0,input,pOrder);
polynom_data = polyval(PP,(min(Pressures):max(Pressures)));

if  max(diff(polynom_data))>=0 &  min(diff(polynom_data))<=0 & sum(abs(diff(input)))>0
    x=diff(polynom_data);
    y=(min(Pressures):max(Pressures));
    jj=0;
    for ii=2:length(polynom_data)-1
        if x(ii-1)<=0 & x(ii)>=0
            jj=jj+1;
            icp(jj)=interp1( [x(ii-1) x(ii)+eps],[y(ii-1) y(ii)],0,'linear') - pBias;
        end
    end
    
    
    if jj==1 
        ICP = mean(icp);
        Factor = 1;
    end
    
    if jj > 1
        ICP = mean(icp);
        Factor = 2;
    end
else
    ICP = 0;
    Factor = 0; 
end

polynom_data = polyval(PP,(0:max(Pressures)));
pFit = polyval(PP, Pressures0);
cc = corrcoef(pFit, input);
Q = cc(2)^2;

