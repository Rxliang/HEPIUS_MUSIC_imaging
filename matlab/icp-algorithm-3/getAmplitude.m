function Amp = getAmplitude(data, fs)

% This function uses the fact that if y = Asin(x), then sqrt(2)*std(y) = A
% For a given segment of length data it calculates the best guess of A
% divides data into overlapping windows and takes the median A of the
% windows

L = length(data);
% we want a minimum of 5 cycles in data - average cycle length is 3-4s.
winLen = round(L/3); %4*fs;
b = 1;
for a = 1:round(winLen/4):L-round(winLen*3/4)-1
    dataWin = detrend(data(a:a+winLen-1));
    dataAmp(b) = sqrt(2)*std(dataWin);
    b = b + 1;
end
if b == 1
    dataAmp(b) = sqrt(2)*std(dataWin);
end
Amp = median(dataAmp);