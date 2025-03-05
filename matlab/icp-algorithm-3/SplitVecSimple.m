function splits = SplitVecSimple(dataIn)

d = diff(dataIn);
f = find(d > 1);
if isempty(f)
    splits = {dataIn};
else
    splits = {zeros(length(f) + 1,1)};
    startSplit = 1;
    f = [f; length(dataIn)];
    for a = 1:length(f)
        splits{a,:} = dataIn(startSplit:f(a));        
        startSplit = f(a) + 1;
    end
end
