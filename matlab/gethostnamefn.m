function hn = gethostnamefn
    
  if ispc | isunix
    [dum, hnPre] = unix('hostname');
  end
  
  hn = strtrim(hnPre);
  
end

