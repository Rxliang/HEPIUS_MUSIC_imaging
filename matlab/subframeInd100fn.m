function iiCell = subframeInd100(sgp,prf)

%extracted from spectralDoppler

%subframeInd100
%partition the data chunk into subrames, one for each line:
%splits up data indices in a frame interval
%R = interpolation factor
R=sgp.R;
TFrame = sgp.TFrame;


ii = round(linspace(1/R,1,R)'*prf*TFrame);
II = [[1;ii(1:end-1)+1],ii];
for k=1:R
    iiCell{k} = II(k,1):II(k,2);
end
end