function [sz, exst] = filesize(filename);

filename

fid = fopen(filename,'r');
if fid < 0
  sz = 0;
  exst = 0;
  return
else
  exst = 1;
end

if feof(fid)
  sz=0;
else
  fseek(fid, 0, 'eof');
  sz = ftell(fid);
  fclose(fid);
end
