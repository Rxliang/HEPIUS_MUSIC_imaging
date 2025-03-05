function [rem, ext] = getfileminusext(filename)

% removes file extension

ind = findstr(filename, '.');

if ~isempty(ind)
  ind = ind(end);
  rem = filename(1:ind-1);
  ext = filename(ind+1:end);
else
  rem=filename;
end

end

