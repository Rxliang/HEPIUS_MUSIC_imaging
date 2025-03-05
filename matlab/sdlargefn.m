%function varargout = sdlargefn(varargin)
function cdiIm = sdlargefn(varargin)

persistent fid
persistent cnt
persistent tm
persistent fileCount
persistent xAx_mm
persistent zAx_mm
persistent fifo

%extclockfn;

PData = evalin('base', 'PData');
PRF = evalin('base', 'evParms.ev.acqPRF');
dopPRF = evalin('base', 'evParms.ev.dopPRF');
nDopPRIUsed = evalin('base', 'evParms.ev.nDopPRIsUsed');
ttna = evalin('base', 'evParms.ev.ttna_microsec')*1e-6;
freq = evalin('base', 'evParms.Trans.frequency');
iqFilePre = evalin('base', 'evParms.file.iqFileSDLargePre');
framesPerFile = evalin('base', 'evParms.file.SDLargeFramesPerFile');
PDataThis = evalin('base', 'PData(evParms.gate.PDataIndSDLarge)');
largeSDSave = evalin('base', 'evParms.largeSDParms.largeSDSave');

if isempty(fileCount)
  fileCount=1;
end
    
IQdata = squeeze(varargin{1});
IQdata = IQdata(:,:,1:nDopPRIUsed);
sz = size(IQdata);

if isempty(cnt)
  fifoLen = 20; % frames
  fifo = dsp.AsyncBuffer(fifoLen);

  cnt=1;
  tm = zeros(framesPerFile,1);
  disp('sdlargefn init');
  dopPRF = -dopPRF; % init signal to mex file
  
  gate = evalin('base', 'evParms.gate');
  origin_mm = gate.lambda_mm*PDataThis.Origin;
  xDelta_mm = gate.lambda_mm*PDataThis.PDelta(1);
  zDelta_mm = gate.lambda_mm*PDataThis.PDelta(3);  
  
  xAx_mm = origin_mm(1):xDelta_mm:origin_mm(1)+(xDelta_mm*(sz(2)-1));
  zAx_mm = origin_mm(3):zDelta_mm:origin_mm(3)+(zDelta_mm*(sz(1)-1));
end

% prepare CDI overlay

np=50;
%np=100;
%np=120;
%np = nDopPRIUsed;

iqmThis = reshape(IQdata(:,:,1:np), [sz(1:2) 1 np]); 

cdiIm = abs(cdiprocoffline(iqmThis));

mx = maxall(cdiIm);
write(fifo, mx);

fifoVals = read(fifo);
medianMax = median(fifoVals);
  
write(fifo, fifoVals);

cdiIm = cdiIm/medianMax*255;

%varargout = cdiIm;

if 0 % plot within own window
  
if cnt == 1
figure(10)
end

imagesc(xAx_mm, zAx_mm,cdiIm);

if cnt==1
  colormap(gray)
  xlabel('mm')
  ylabel('mm')
end

drawnow

end

if largeSDSave  

  iqFileRe = [iqFilePre '_' num2str(fileCount) '.re'];
  iqFileIm = [iqFilePre '_' num2str(fileCount) '.im'];

  dateNum = datenum(now);

  if cnt==framesPerFile
    closeFilesAfterWrite=1;
  else
    closeFilesAfterWrite=0;
  end
  tic;
  saveiq(IQdata, PRF, dopPRF, ttna, nDopPRIUsed, iqFileRe, iqFileIm, dateNum, ...
         PDataThis.Origin, PDataThis.PDelta, freq, closeFilesAfterWrite);
  tm(cnt)=toc;
  
  cnt=cnt+1;

  if cnt==framesPerFile+1
    cnt=1;
    fileCount = fileCount+1;
  end

  assignin('base', 'tm', tm);
else
  cnt=cnt+1;    
end




%if isempty(fid)

   
%  dateStr = datestr(now, 'yyyymmdd_HHMM');
%  mfile = mfilename;
%  vtData = [getenv('VITTAMED_DATA') '/nicp/'];
%  fid = fopen([vtData mfile '_' dateStr], 'wb');
%  fwrite(fid, size(IQdata), 'uint16');
%  fwrite(fid, prf, 'double');  
%end

%tic
%parfor i = 1:1 
%  fwrite(fid, IQdata, 'double');
%end
%toc

%size(IQdata)

