function UI = inithvfn(UI, u, VInit)
   
for profNo = 1:length(VInit)
  UI(u).Statement = ['[result,hv] = setTpcProfileHighVoltage(' ...
                     num2str(VInit(profNo)) ',' num2str(profNo) ');'];
  u=u+1;
  sldrStr = ['''hv' num2str(profNo) 'Sldr'''];
  UI(u).Statement = ['hvSldr = findobj(''tag'',' sldrStr ');' ...
                     'set(hvSldr,''Value'', ' num2str(VInit(profNo)) ...
                     ');'];
  u=u+1;
  valStr = ['''hv' num2str(profNo) 'Value'''];
  UI(u).Statement = ['hvValue = findobj(''tag'', ' valStr '); '...
                     'set(hvValue,''String'', ' ...
                     num2str(VInit(profNo), '%.1f') ');'];
  u=u+1;
end