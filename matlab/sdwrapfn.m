function varargout = sdwrapfn(varargin)

persistent instanceNumber enabledByDefault

if ~isempty(instanceNumber) & enabledByDefault
  varargout = feval(['sdfn_' num2str(instanceNumber)], varargin);
  return
end
evParms = evalin('base', 'evParms');

if isempty(instanceNumber)
  mfile = mfilename;
  instanceNumber = str2num(mfile(end));
  if ~isfield(evParms.SD, 'enableInstances');
    enabledByDefault=1;
  else
    enabledByDefault=0;
  end
end

if enabledByDefault | ismember(instanceNumber, evParms.SD.enableInstances) 
  % note, old Vantage versions passed interbuffer as complex single 
  % argument. later v4s pass as real and imag arg 1 and 2!
    %  varargout{1} = feval(['sdfn_' num2str(instanceNumber)], varargin{1});
  
  IQin = complex(varargin{1}, varargin{2});
  varargout{1} = feval(['sdfn_' num2str(instanceNumber)], IQin);
else
  varargout{1} = {};
end



