function done = pendingcontrolprocfn(sleepTime_s)

Control = evalin('base', 'Control');

%length(Control)
%isempty(Control(1).Command)
%isempty(Control(1).Parameters)

done = (length(Control)==1 & isempty(Control(1).Command) & ...
           isempty(Control(1).Parameters));



