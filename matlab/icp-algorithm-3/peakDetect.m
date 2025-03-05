function numPeaks = peakDetect(data, pctChange)

% Function to perform large peak detection 
L = length(data);

% peak detection
firstDiff = diff(data);
c = firstDiff(1:end-1).*firstDiff(2:end);
f = find(c <= 0) + 1;  % zero crossings
if f(1) == 1
    f(1) = [];
end
d = find( firstDiff(f-1) >= 0 & firstDiff(f) < 0);
dataPeaks = f(d);

% if isempty(dataPeaks)
%     numPeaks = 0;
%     return
% end

% find troughs
d = find( firstDiff(f-1) < 0 & firstDiff(f) >= 0);
dataTroughs = f(d);

% line up data
if dataPeaks(1) < dataTroughs(1)
    dataPeaks(1) = [];
end

minLen = min(length(dataTroughs), length(dataPeaks));
amp = data(dataPeaks(1:minLen)) - data(dataTroughs(1:minLen));
amp(amp == 0) = [];

thresh = prctile(amp,95)*pctChange;

numPeaks = length(find(amp > thresh));

