function h = linept2pt(p1v,p2v)

if min(size(p1v)) == 1
  p1v = p1v(:).';
end

if min(size(p2v)) == 1
  p2v = p2v(:).';
end


for i = 1:size(p1v,1)
  p1 = p1v(i,:);
  p2 = p2v(i,:);  

  if length(p1)==3
    h(i) = line([p1(1) p2(1)], [p1(2) p2(2)], [p1(3) p2(3)]);
  else
    h(i) = line([p1(1) p2(1)], [p1(2) p2(2)]);  
  end

end
