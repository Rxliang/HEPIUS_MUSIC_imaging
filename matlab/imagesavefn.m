function imagesavefn(imIn)

persistent cnt 
persistent outFile


maxFrames = 500;
instanceNumber = 1; % B

return

if ~exist('cnt', 'var') | isempty(cnt)
  cnt=0;
  assignin('base', 'imIn', imIn);
  evalin('base', ...
  ['imSet.im = zeros(size(imIn,1), size(imIn,2), ' num2str(maxFrames) ...
   ', ''like'', imIn);']);
  
  datenumStart = datenum(now);
  dateStr = datestr(datenumStart, 'yyyymmdd_HHMM');
  outFilePre = evalin('base', 'outFilePre');
  outFile = [outFilePre '_' dateStr];
  PData = evalin('base', 'PData');
  spec.PData = PData(1);
  spec.lambda_mm = evalin('base', 'lambda_mm');
  firstRun = 1;
  h5File = writebmodeh5fxffn(spec, imIn, 'B', outFile, 'single', ...
                             instanceNumber, firstRun);
end

cnt=cnt+1;
dnow = datenum(now);
dispdiv('cnt', cnt, 500);
if 0
  ind = mod(cnt-1, maxFrames)+1
  assignin('base', 'imIn', imIn);
  assignin('base', 'ind', ind);
  evalin('base', ['imSet.im(:,:,ind) = imIn;']);
  
  evalin('base', ['imSet.datenum(ind) = ' num2str(dnow) ';']);
end

spec.datenumFrame = dnow;
tic
size(imIn)
h5File = writebmodeh5fxffn(spec, imIn, 'B', outFile, 'single', ...
                           instanceNumber, 0);
toc
end

