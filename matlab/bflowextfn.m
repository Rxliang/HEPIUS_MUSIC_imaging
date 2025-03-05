function outBF = bflowextfn(imIn)

persistent cnt 
persistent outFile
persistent outBFPre
persistent h

maxFrames = 500;
instanceNumber = 1; % B

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
  %h5File = writebmodeh5fxffn(spec, imIn, 'B', outFile, 'single', ...
  %                           instanceNumber, firstRun);
  outBFPre = [];
  figure(2)
  clf
  h = imagesc(imIn);
  colormap(gray)
end

cnt=cnt+1;
frameIndex = mod(cnt-1,3)+1;

outBFPre(:,:,frameIndex) = imIn;

if frameIndex == 3
  outBF = mean(outBFPre(:)) * 0.5* ((outBFPre(:,:,3)-outBFPre(:,:,2)) + (outBFPre(:,:,2)-outBFPre(:,:,1)));
  set(h, 'CData', outBF);
else
  outBF = zeros(size(imIn));
end
     
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
%tic
%size(imIn)
%h5File = writebmodeh5fxffn(spec, imIn, 'B', outFile, 'single', ...
%                           instanceNumber, 0);
                           %toc
end

