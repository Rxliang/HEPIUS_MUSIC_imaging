% Script to analyze (compare) data in the studyStore data structure

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

% % User selects which files to compare
numFiles = length(studyStore{1, folderSelection}.data);

[intextSelection1, ok] = listdlg('PromptString','Select in or ext for 1st file:', ...
    'SelectionMode','single','ListString', {'Internal', 'External'});

if ~ok
    return
end

count = 1;
fileIndex = [];
for a = 1:numFiles
    if ~isempty(studyStore{1, folderSelection}.data(1, a).int)
        fileString{count} = num2str(a);
        fileIndex(count) = a;
        count = count + 1;
    end
end

[fileSelection1, ok] = listdlg('PromptString','Select 1st file:', ...
    'SelectionMode','single','ListString', fileString);

if ~ok
    return
end

[intextSelection2, ok] = listdlg('PromptString','Select in or ext for 2nd file:', ...
    'SelectionMode','single','ListString', {'Internal', 'External'});

if ~ok
    return
end

[fileSelection2, ok] = listdlg('PromptString','Select 2nd file:', ...
    'SelectionMode','single','ListString', fileString);

if ~ok
    return
end

h(1) = subplot(2,1,1);
f1 = fileIndex(fileSelection1);
if intextSelection1 == 1
    imagesc(log10(studyStore{1, folderSelection}.data(1, f1).int.averagedSono));
else
    imagesc(log10(studyStore{1, folderSelection}.data(1, f1).ext.averagedSono));
end

h(2) = subplot(2,1,2);
f2 = fileIndex(fileSelection2);
if intextSelection1 == 1
    imagesc(log10(studyStore{1, folderSelection}.data(1, f2).int.averagedSono));
else
    imagesc(log10(studyStore{1, folderSelection}.data(1, f2).ext.averagedSono));
end

linkaxes(h, 'x')