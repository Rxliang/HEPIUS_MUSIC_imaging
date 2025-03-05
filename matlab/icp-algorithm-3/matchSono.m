function [KLD] = matchSono(im1, im2, numMerge)
%UNTITLED2 [match] = matchSono(im1, im2)
%   Matches vertical histograms using K-L divergence

%bins = 1:.25:20;
% binMin = min(min(min(im1)), min(min(im2)));
% binMax = max(max(max(im1)), max(max(im2)));
% bins = linspace(binMin, binMax, 20);
bins = linspace(0, 1, 20);

[numRows, numCols] = size(im1);

for a = 1:numCols-numMerge+1
    data1 = im1(:,a:a+numMerge-1);
    data2 = im2(:,a:a+numMerge-1);
    
    p1 = hist(data1(:), bins);
    p1 = p1'./sum(p1) + eps;
    
    p2 = hist(data2(:), bins);
    p2 = p2'./sum(p2) + eps;
    
    KLD(a) = kldiv(bins', p1, p2);
end


end

