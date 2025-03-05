function Callback = setprocparmfn(parmName, procIndName, updateVar)

evParms = evalin('base', 'evParms');

procStr = ['procIndThis = evParms.proc.' procIndName ';'];
procIndNameStr = ['procIndName = ''' procIndName ''';'];
parmStr = ['parmName = ''' parmName ''';'];
%updateStr = [updateVar ' = num2str(UIValue);']; 
updateStr = [updateVar ' = UIValue;']; 
Callback = text2cell('%macro1');
Callback{2} = procStr;
Callback{3} = procIndNameStr;
Callback{4} = parmStr;
Callback{5} = updateStr;
Callback = [Callback text2cell('%macro2')];

return

%macro1
evParms = evalin('base', 'evParms');
% need space here
%macro1

%macro2
Process = evalin('base', 'Process');

for k = 1:2:length(Process(procIndThis).Parameters) if strcmp(Process(procIndThis).Parameters{k}, parmName) Process(procIndThis).Parameters{k+1} = UIValue; end; end

assignin('base','Process',Process);
Control = evalin('base','Control');
Control.Command = 'set&Run';
Control.Parameters = {'Process', procIndThis, parmName , UIValue};
assignin('base','Control', Control);
assignin('base','evParms', evParms);
disp(['Set ' parmName ' of process ' procIndName ' (' num2str(procIndThis) ') to ' num2str(UIValue)]);
%macro2