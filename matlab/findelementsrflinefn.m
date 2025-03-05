function evParms = findelementsrflinefn(evParms, Trans)

evParms.BLineWT.elements = [];
evParms.BLineWT.channels = [];

for q = 1:evParms.BLineWT.numOrigin
  [dum, indCenterElement] = findclosestinvec(Trans.ElementPos(:,1), ...
                                             evParms.gate.x(q));
  evParms.BLineWT.elements{q} =  ceil(indCenterElement- ...
                                   evParms.BLineWT.numChannels/2):floor(indCenterElement+ ...
                                                    evParms.BLineWT.numChannels/2);
  evParms.BLineWT.elements{q} = unique(min(evParms.BLineWT.elements{q}, ...
                                        evParms.ev.numElementsUsed));
  evParms.BLineWT.elements{q} = unique(max(evParms.BLineWT.elements{q}, 1)); ...

  if isfield(Trans, 'ConnectorES')
      evParms.BLineWT.channels{q} = Trans.ConnectorES(evParms.BLineWT.elements{q});
  else
      evParms.BLineWT.channels{q} = Trans.Connector(evParms.BLineWT.elements{q});
  end
  
end


end

