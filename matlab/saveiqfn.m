function saveiqfn(IQdata, evParms, PData, iqFileRe, iqFileIm, dateNum, ...
    instanceNumberIn, saveImFlag)

persistent cnt
persistent instanceNumber

if nargin < 8
  saveImFlag = 0;
end

if nargin < 7 & isempty(instanceNumber)
  mfile = mfilename;
  instanceNumber = str2num(mfile(end));
else
  instanceNumber = instanceNumberIn;    
end

if isempty(cnt)
  cnt=1;
 % dopPRF = -evParms.ev.dopPRF; % init signal to mex file
  dopPRF = evParms.ev.dopPRF; % init signal no longer needed
else
  dopPRF = evParms.ev.dopPRF;
end
ttna = evParms.ev.ttna_microsec*1e-6;
PRF = evParms.ev.acqPRF;

%IQdata = squeeze(varargin{1});

% new file for each invocation
closeFileAfterWriting=1;

freq = evParms.Trans.frequency;

calledFn = ['saveiq_' num2str(instanceNumber)];

%disp(['saveiqfn: in instance:' mfilename]);

tic
fileVersion=3;
feval(calledFn, single(IQdata), PRF, dopPRF, ttna, evParms.ev.nDopPRIs, ...
       iqFileRe, iqFileIm, dateNum, PData.Origin, PData.PDelta, freq, ...
       closeFileAfterWriting, instanceNumber, char(evParms.this.mode), ...
       double(evParms.this.TWFreq_MHz), ...
       uint16(evParms.gate.lateralSubSampleIndex{instanceNumber}), ...
       uint16(evParms.gate.rangeSubSampleIndex{instanceNumber}), ...
       double(fileVersion));

if saveImFlag
    imFile = [getfileminusext(iqFileRe) '.img'];
    instanceNumberUsed=1;
    calledFn = 'saveim';
    imFileVersion=1;
    PDataBase = evalin('base', 'PData');                    
    PDataThis = PDataBase(evParms.gate.PDataIndB);
    pdeltaIm = evalin('base','Resource.DisplayWindow(1).pdelta');  
    hIm = evalin('base','Resource.DisplayWindow(1).imageHandle');    
    imBuf = get(hIm, 'CData');        
    PDataThis.PDelta(2) = pdeltaIm; % steal this y-storage space
    feval(calledFn, single(imBuf), PRF, dopPRF, ttna, ...
      evParms.ev.nDopPRIs, ...
      imFile,  dateNum(1), PDataThis.Origin, ...
      PDataThis.PDelta, freq, ...
      closeFileAfterWriting, instanceNumberUsed, ...
      char(evParms.this.mode), ...
      double(evParms.this.TWFreq_MHz), ...
      double(imFileVersion));


end



toc;

cnt = cnt + 1;
end

