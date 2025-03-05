% script to cycle through images

if ~exist('studyStore', 'var')
    errordlg('Please load in studyStore data', 'Error');
    return
end

% User selects which study folder to compare
numFolders = length(studyStore);
for a = 1:numFolders
    listString{a} = studyStore{a}.name;
end

[folderSelection, ok] = listdlg('PromptString','Please select a study folder:', ...
    'SelectionMode','single','ListString',listString);

if ~ok
    return
end

numFiles = length(studyStore{1, folderSelection}.data);

count = 1;
fileIndex = [];
for a = 1:numFiles
    if ~isempty(studyStore{1, folderSelection}.data(1, a).int)
        fileString{count} = num2str(a);
        fileIndex(count) = a;
        count = count + 1;
    end
end

bPlot = 1;
match1 = [];

sonoObj = sono;

for a = 1:count-1
    im1 = studyStore{1, folderSelection}.data(1, fileIndex(a)).int.averagedSono;
    im1Full = studyStore{1, folderSelection}.data(1, fileIndex(a)).int.AR;
    env1 = studyStore{1, folderSelection}.data(1, fileIndex(a)).int.env;
%    im1 = im1./nansum(nansum(im1));
    im1 = im1./max(max(im1));
    im1Full = im1Full./max(max(im1Full));
    if bPlot
        subplot(1,2,1)
        imagesc(log10(imStretch(im1)));
        title(['Internal: ADC ' num2str(fileIndex(a))]);
    end
    
%     [nRows, nCols] = size(im1);
%     mom1 = zeros(1, nCols);
%     for b = 1:nCols
%         aa = isnan(im1(:,b));
%         yVal = find(~aa, 1, 'first');
%         mom1(b) = trapz( (1:nRows-yVal+1)'.^1.*im1(yVal:end,b) );
%     end
%     m1(a) = sum(mom1);
    
    
    im2 = studyStore{1, folderSelection}.data(1, fileIndex(a)).ext.averagedSono;
    im2Full = studyStore{1, folderSelection}.data(1, fileIndex(a)).ext.AR;
    env2 = studyStore{1, folderSelection}.data(1, fileIndex(a)).ext.env;
  %  im2 = im2./nansum(nansum(im2));
    im2Full = im2Full./max(max(im2Full));
    
    if bPlot
        subplot(1,2,2)
        imagesc(log10(imStretch(im2)));
        title(['External: ADC ' num2str(fileIndex(a))]);
    end
    
%     [nRows, nCols] = size(im2);
%     mom2 = zeros(1, nCols);
%     for b = 1:nCols
%         aa = isnan(im2(:,b));
%         yVal = find(~aa, 1, 'first');
%         mom2(b) = trapz( (1:nRows-yVal+1)'.^1.*im2(yVal:end,b) );
%     end
%     m2(a) = sum(mom2);
    
    %kld = matchSono((im1), (im2), 2);
    %match1(a) = mean(kld);
    %disp(mean(kld))
    %mom1 = sonoObj.spectralMoment(flipud(im1Full), size(im1Full,1) - env1, 1);
    %mom2 = sonoObj.spectralMoment(flipud(im2Full), size(im2Full,1) - env2, 1);
    mom1 = sonoObj.spectralMoment(im1Full, env1, 1);
    mom2 = sonoObj.spectralMoment(im2Full, env2, 1);
    matches = sonoObj.matchWaveforms(mom1, mom2, env1, env2, 1);
    disp(['Num Good Cycles: ' num2str(matches.numGoodPeriods)]);
    disp(['SSE: ' num2str(matches.SSE)]);

    if bPlot
        pause;
    end
end


