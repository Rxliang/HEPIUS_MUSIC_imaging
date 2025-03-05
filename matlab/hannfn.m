function w = hannfn(M)

% ccrma.stanford.edu/~jos/sasp/Matlab_Hann_Window.html

w = .5*(1 - cos(2*pi*(0:M-1)'/(M-1)));