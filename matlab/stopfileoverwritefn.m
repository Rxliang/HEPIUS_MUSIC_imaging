function outFile = stopfileoverwritefn(outFilePre, fileExt)
    
% prevent overwrites
fileFixStr = '';
cntFix=0;
cont=1;

while cont

  cond = 0;
  for q = 1:length(fileExt);
     imFile = [outFilePre fileFixStr fileExt{q}];
     cond = cond | exist(imFile, 'file');
  end
    
  if cond
      disp(['*** Warning: ' imFile ' exists!']);
      cntFix=cntFix+1;
      fileFixStr = ['_fix' num2str(cntFix)]; 
  else
      cont=0;
  end            
end

outFile = [outFilePre fileFixStr];

end




