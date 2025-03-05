function lslrt(direPre)

if isunix
  v = ver('matlab');
%  dire = escapespace(direPre);
  dire = direPre;
  ind = findstr(v.Version, '.');
  vMain = str2num(v.Version(1:ind-1));
  vSub = str2num(v.Version(ind+1:end));
  if (vMain == 7 & vSub >= 7) | (vMain > 7)
    ls('-lrt', [dire]);
  else
    ls('-lrt', ['"' dire '"']);
  end
else
  fix1 = replacechar(direPre,'/','\');
  evalStr = ['! dir "' fix1 '"'];
  eval(evalStr);
end

