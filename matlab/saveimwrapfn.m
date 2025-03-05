function saveimwrapfn(imBuf)

global semaphoreKey; semaphoreKey = 24537;
global recordDataFileName; 
global recordData;

instanceNumber = 1;
evParms = evalin('base', 'evParms');
semaphoreOpenedState = evParms.state.semaphoreOpenedState;
PData = evalin('base', 'PData');                    
PDataThis = ...
    PData( evParms.gate.PDataIndSDLarge);

imLen = PDataThis.Size(1)*PDataThis.Size(2);

im = reshape(imBuf(1:imLen), PDataThis.Size(1), PDataThis.Size(2));

if isempty(instanceNumber)
  mfile = mfilename;
  instanceNumber = str2num(mfile(end));
end

NframesIn = []; % read from evParms

[dataFileName]=saveimwithtimestampfn(im, PDataThis, ...
                                     NframesIn, instanceNumber,...
                                     semaphoreOpenedState);
 
