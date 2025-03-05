function [val, ind] = findclosestinvec(vec, inVals)

if length(inVals) == 1
  absDiff =  abs(vec-inVals);
  [dum, ind] = min(absDiff);
  val = vec(ind);
else
  inVals = inVals(:)';
  vec = vec(:);
  inValsRep  = repmat(inVals, length(vec),1);
  vecRep  = repmat(vec, 1, length(inVals));
  absDiff =  abs(vecRep-inValsRep);
  [dum, ind] = min(absDiff);
  val = reshape(vec(ind,:),size(inVals));
end