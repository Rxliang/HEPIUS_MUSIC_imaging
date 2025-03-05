function central = findcentraln(n, inVec)

N = length(inVec);

midInd = floor(N/2)+1;

rLen = round(n/2);
fLen = floor(n/2);

if ~iseven(N)
  if ~iseven(n)
    central = inVec(midInd-fLen:midInd+fLen);
  else
    central = inVec(midInd-fLen+1:midInd+fLen);
  end
else
  if ~iseven(n)
    central = inVec(midInd-fLen:midInd+fLen);
  else
    central = inVec(midInd-fLen:midInd+fLen-1);
  end
  
end








