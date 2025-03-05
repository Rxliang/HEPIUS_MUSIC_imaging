function bool = dispdiv(counter, cnt, interval)
q = single(cnt) / interval;

if q == round(q)
  disp([counter ' = ' num2str(cnt)]);
  bool = 1;
else
  bool = 0;
end
