
inFilePre = '20240125_195309_s1_t1_sg3_i1';
inFilePre = '20240125_230210_s1_t1_sg1_i1';

matInPath = 'C:\Users\223022005\Documents\music\data\iqdata\semaphoretest\';

matInPath = 'C:\Users\Administrator\Documents\music\data\string';
inFilePre = 'string20x1020240330_130432_sdl_s4_sg1_fix1'
inFilePre = 'string20x1020240330_133818_sdl_s5_sg3'

fileVersion = 1;
overWrite = 1;
[outFile, iSet] = readimfn(inFilePre, matInPath, fileVersion, ...
                                     overWrite)

figure(10)
clf
for q = 1:size(iSet.im,3);
    imagesc(iSet.xAx_mm, iSet.zAx_mm, iSet.im(:,:,q))
    drawnow
end

