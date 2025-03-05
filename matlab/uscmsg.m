function uscmsg(execState)

if strcmp(execState, 'production')
  disp('*** running unsupported case in production run *** ');
  pausede
end
