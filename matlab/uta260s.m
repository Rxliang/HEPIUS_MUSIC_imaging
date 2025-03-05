function uta260sOut = uta260s

vtHome = getenv('MUSIC_HOME');
inFile = [vtHome '/docs/uta260s.xlsx'];
lslrt(inFile)
tab = readtable(inFile, 'ReadVariableNames', 0);
%tab = xlsread(inFile);
cell = table2cell(tab);

numBanks = size(cell,2)-1;
numRows =  size(cell,1)/2;

cnt = 0;
uta260sOut.pinNumber = {};
uta260sOut.pinName = [];

endChar = char(160);

for i = 1:numBanks
  for j = 1:numRows
    cnt=cnt+1;
    rowInd = 2*(j-1)+1;
    uta260sOut.pinNumber{cnt} = replacechar(cell{rowInd,i+1}, endChar,''); 
    uta260sOut.pinName{cnt} = replacechar(cell{rowInd+1,i+1}, endChar,''); 
  end
end


