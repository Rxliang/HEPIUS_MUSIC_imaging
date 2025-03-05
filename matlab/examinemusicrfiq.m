if ~exist('rfiq', 'var') | isempty(rfiq)
  rfiq = load('../data/20231109_Plane_noavg_arr3.mat');
end

figure(1)
clf
imagesc(log(rfiq.ImgData{1}(:,:,1,1)))
colormap gray

figure(2)
clf
imagesc(log(rfiq.ImgData{1}(:,:,1,2)))
colormap gray

d1 = rfiq.RcvData{1};
i1 = d1(:,:,1);
% see which channels are active
s1i1 = sum(i1,1);
indCh = find(s1i1)

i1Sel = i1(:,indCh);

figure(3)
plot(i1Sel)


d2 = rfiq.RcvData{2};
i2 = d2(:,:,1);
% see which channels are active

i2Sel = i2(:,indCh);

figure(4)
plot(i2Sel)


d3 = rfiq.RcvData{3};
i3 = d3(:,:,1);
% see which channels are active

i3Sel = i3(:,indCh);

figure(5)
plot(i3Sel)


d4 = rfiq.RcvData{4};
i4 = d4(:,:,1);
% see which channels are active

i4Sel = i4(:,indCh);

figure(6)
plot(i4Sel)
