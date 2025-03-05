function outStr = replacechar(inStr, ch, chNew)

if length(chNew)==1
  ind = find(inStr==ch);  
  outStr=inStr;
  outStr(ind) = chNew;
else
  ind = findstr(inStr,ch);
  if isempty(ind)
    outStr = inStr;
    return
  end
  
  outStr = [];
  cnt=1;
  for i = 1:length(ind)
    outStr = [outStr inStr(cnt:ind(i)-1) chNew];
    cnt = ind(i)+length(chNew)-1;
  end
  outStr = [outStr inStr(ind(end)+1:end)];
  
end
