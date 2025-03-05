function axc = axbound2centre(axb)

axc = axb(1:end-1)+diff(axb)/2;

