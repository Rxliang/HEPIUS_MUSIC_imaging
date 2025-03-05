% run readmmodever1 1st

figure(1)
clf

numLines = size(MSet.M,2);

for q = 149:numLines
  plot(MSet.zAx_mm, MSet.M(:,q).');
  pausede
end
