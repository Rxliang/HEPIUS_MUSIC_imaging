function out = titstrfn(in, len)
  if nargin < 2
    len = 80;
  end
  in = replacechar(in, '_', '\_');

  
  cnt = 1;
  out = [];
  ln=1;
  L = length(in);
  
  if len > L
    len=L;
  end
  
  while cnt < L
    en = min(L, cnt+len-1);
    seg = in(cnt:en);
    if en < L 
      spInd = find((seg==' ' | seg=='.' | seg == '-'),1,'last');
      if ~isempty(spInd)
          seg= seg(1:spInd);
      end
    end
    
    cnt = cnt+length(seg);
       
    out{ln} = seg;
    ln=ln+1;
  end
    
end

