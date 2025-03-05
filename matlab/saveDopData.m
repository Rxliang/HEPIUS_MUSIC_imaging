function saveDopData(fname)

Event = evalin('base','Event');
RcvData = evalin('base','RcvData');
IData = evalin('base','IData');
QData = evalin('base','QData');
ImgData = evalin('base','ImgData');
Receive = evalin('base','Receive');
Trans = evalin('base','Trans');
TX = evalin('base','TX');
ReconInfo = evalin('base','ReconInfo');
RcvProfile = evalin('base','RcvProfile');
PData = evalin('base','PData');
dop = evalin('base','dop');
P = evalin('base','P');
Recon = evalin('base','Recon');
% cdi = evalin('base','cdi');

% savefast(fname,'Event','IData','QData','Receive','Trans','TX','PData','RcvProfile','ReconInfo','ImgData','RcvData','dop','P','cdi')
savefast(fname,'Event','IData','QData','Receive','Trans','TX','PData','RcvProfile','ReconInfo','ImgData','RcvData','dop','P','Recon', 'Resource', 'SeqControl')

% JSM added Resource, SeqControl so EventAnalysisTool might work


end

