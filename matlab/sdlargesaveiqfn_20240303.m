function sdlargesaveiqfn(IQdataRe, IQdataIm)

persistent frameCount prevOutFilePrefix dateNumVec frameStore

IQdata = IQdataRe + j*IQdataIm;

% sz 406   224     1   226
keyboard

if ~exist('frameStore', 'var') | isempty(frameStore)
  frameStore = zeros(346, 187, 226, 100, 'single');
end

if ~exist('frameCount', 'var') | isempty(frameCount)
  frameCount = 0;
end

evParms = evalin('base', 'evParms');

if ~isfield(evParms.state, 'recording') | ...
      isempty(find(evParms.state.recording))

  if frameCount > 0
    PData = evalin('base', 'PData');                    
    PDataThis = PData(evParms.gate.PDataIndSDLarge);

    currentOutputFilePrefix = evalin('base', ...
                                 'currentOutputFilePrefix');
    iqFileRe = [currentOutputFilePrefix '_sdl.re'];
    iqFileIm = [currentOutputFilePrefix '_sdl.im'];

    saveiqfn(frameStore, evParms, PData, iqFileRe, iqFileIm, dateNumVec);

    frameCount = 0;
    dateNumVec = [];
%    frameStore = [];

  end
  return
end

frameCount = frameCount+1
dateNumVec(frameCount) = datenum(now);
size(IQdata)
frameStore(:, :, :, frameCount) = squeeze(IQdata);


%if ~exist('prevOutFilePrefix', 'var')
%  prevOutFilePrefix = currentOutputFilePrefix;
%end



end

