% Batchscript to go through Aarau data and extract Vd for Int (2) and Ext
% (1)

load('vdStats.mat');
% re-organize to map vdStats analysis for TCD3D Vd analysis
numFolders = size(vdStatsInt, 2);
for a = 1:numFolders
    vdStats.v2group{a} = vdStatsInt{a}.vd;
    vdStats.v1group{a} = vdStatsExt{a}.vd;
    vdStats.presGroup{a} = unique(vdStatsInt{a}.pressure);
end

% organize by pressure
err1 = NaN;
for a = 1:numFolders
    tmp = vdStats.v1group{a}(find(vdStats.presGroup{a} == 0));
    if ~isempty(tmp), vdStats.v1_0(a) = tmp; else vdStats.v1_0(a) = err1; end
    tmp = vdStats.v1group{a}(find( (vdStats.presGroup{a} == 4)|(vdStats.presGroup{a} == 3)));
    if ~isempty(tmp), vdStats.v1_4(a) = tmp; else vdStats.v1_4(a) = err1; end
    tmp = vdStats.v1group{a}(find( (vdStats.presGroup{a} == 8)|(vdStats.presGroup{a} == 9)));
    if ~isempty(tmp), vdStats.v1_8(a) = tmp; else vdStats.v1_8(a) = err1; end
    tmp = vdStats.v1group{a}(find(vdStats.presGroup{a} == 12));
    if ~isempty(tmp), vdStats.v1_12(a) = tmp; else vdStats.v1_12(a) = err1; end
    tmp = vdStats.v1group{a}(find( (vdStats.presGroup{a} == 16) | (vdStats.presGroup{a} == 15)));
    if ~isempty(tmp), tmp = tmp(1); end
    if ~isempty(tmp), vdStats.v1_16(a) = tmp; else vdStats.v1_16(a) = err1; end
    tmp = vdStats.v1group{a}(find(vdStats.presGroup{a} == 20));
    if ~isempty(tmp), vdStats.v1_20(a) = tmp; else vdStats.v1_20(a) = err1; end
    tmp = vdStats.v1group{a}(find(vdStats.presGroup{a} == 24));
    if ~isempty(tmp), vdStats.v1_24(a) = tmp; else vdStats.v1_24(a) = err1; end
    tmp = vdStats.v1group{a}(find(vdStats.presGroup{a} == 28));
    if ~isempty(tmp), vdStats.v1_28(a) = tmp; else vdStats.v1_28(a) = err1; end
       
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 0));
    if ~isempty(tmp), vdStats.v2_0(a) = tmp; else vdStats.v2_0(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 4));
    if ~isempty(tmp), vdStats.v2_4(a) = tmp; else vdStats.v2_4(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 8));
    if ~isempty(tmp), vdStats.v2_8(a) = tmp; else vdStats.v2_8(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 12));
    if ~isempty(tmp), vdStats.v2_12(a) = tmp; else vdStats.v2_12(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 16));
    if ~isempty(tmp), vdStats.v2_16(a) = tmp; else vdStats.v2_16(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 20));
    if ~isempty(tmp), vdStats.v2_20(a) = tmp; else vdStats.v2_20(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 24));
    if ~isempty(tmp), vdStats.v2_24(a) = tmp; else vdStats.v2_24(a) = err1; end
    tmp = vdStats.v2group{a}(find(vdStats.presGroup{a} == 28));
    if ~isempty(tmp), vdStats.v2_28(a) = tmp; else vdStats.v2_28(a) = err1; end    
end

% get differences
vdStats.v1_0_4 = (vdStats.v1_0 - vdStats.v1_4)./vdStats.v1_0;
vdStats.v1_0_8 = (vdStats.v1_0 - vdStats.v1_8)./vdStats.v1_0;
vdStats.v1_0_12 = (vdStats.v1_0 - vdStats.v1_12)./vdStats.v1_0;
vdStats.v1_0_16 = (vdStats.v1_0 - vdStats.v1_16)./vdStats.v1_0;
vdStats.v1_0_20 = (vdStats.v1_0 - vdStats.v1_20)./vdStats.v1_0;
vdStats.v1_0_24 = (vdStats.v1_0 - vdStats.v1_24)./vdStats.v1_0;
vdStats.v1_0_28 = (vdStats.v1_0 - vdStats.v1_28)./vdStats.v1_0;

vdStats.v2_0_4 = (vdStats.v2_0 - vdStats.v2_4)./vdStats.v2_0;
vdStats.v2_0_8 = (vdStats.v2_0 - vdStats.v2_8)./vdStats.v2_0;
vdStats.v2_0_12 = (vdStats.v2_0 - vdStats.v2_12)./vdStats.v2_0;
vdStats.v2_0_16 = (vdStats.v2_0 - vdStats.v2_16)./vdStats.v2_0;
vdStats.v2_0_20 = (vdStats.v2_0 - vdStats.v2_20)./vdStats.v2_0;
vdStats.v2_0_24 = (vdStats.v2_0 - vdStats.v2_24)./vdStats.v2_0;
vdStats.v2_0_28 = (vdStats.v2_0 - vdStats.v2_28)./vdStats.v2_0;

% reorganize for anova1
vdStats.v1_all = [vdStats.v1_0'; vdStats.v1_4'; vdStats.v1_8'; vdStats.v1_12'; vdStats.v1_16'; vdStats.v1_20'; vdStats.v1_24'; vdStats.v1_28'];
vdStats.v2_all = [vdStats.v2_0'; vdStats.v2_4'; vdStats.v2_8'; vdStats.v2_12'; vdStats.v2_16'; vdStats.v2_20'; vdStats.v2_24'; vdStats.v2_28'];
vdStats.v1_diff_all = [vdStats.v1_0_4'; vdStats.v1_0_8'; vdStats.v1_0_12'; vdStats.v1_0_16'; vdStats.v1_0_20'; vdStats.v1_0_24'; vdStats.v1_0_28'];
vdStats.v2_diff_all = [vdStats.v2_0_4'; vdStats.v2_0_8'; vdStats.v2_0_12'; vdStats.v2_0_16'; vdStats.v2_0_20'; vdStats.v2_0_24'; vdStats.v2_0_28'];

L1 = ones(numFolders,1);
vdStats.L = [L1*0; L1*4; L1*8; L1*12; L1*16; L1*20; L1*24; L1*28];
vdStats.L_diff = [L1*0; L1*4; L1*8; L1*12; L1*16; L1*20; L1*24];
    
idx = 1;
for a = 1:7
    v1_diff(a) = nanmean(vdStats.v1_diff_all(idx:idx+numFolders-1));
    v2_diff(a) = nanmean(vdStats.v2_diff_all(idx:idx+numFolders-1));
    idx = idx+numFolders;
end
pres = [0; 4; 8; 12; 16; 20; 24];

[p,tbl,stats] = anova1(vdStats.v1_diff_all, vdStats.L_diff);
[c,~,~,gnames] = multcompare(stats);

0 - 4
0 - 8
0 - 12
0 - 16
0 - 20
0 - 24
