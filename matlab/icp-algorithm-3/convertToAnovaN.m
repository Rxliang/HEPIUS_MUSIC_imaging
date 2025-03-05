function balOut = convertToAnovaN(bal1, bal2)

% function to convert balance factor structure for use by anovan
bal = [bal1 bal2];
bal = bal(:);

[numPressures1, numSamples1] = size(bal1);
[numPressures2, numSamples2] = size(bal2);
g1a = 1:1:numPressures1;
g1b = 1:1:numPressures2;

g1a = repmat(g1a, 1, numSamples1);
g1b = repmat(g1b, 1, numSamples2);
g1 = [g1a g1b];

L1 = numSamples1*numPressures1;
L2 = numSamples2*numPressures2;
g2 = ones(1, L1);
g2(L1+1:L1+L2) = ones(1, L2)*2;

balOut.g1 = g1;
balOut.g2 = g2;
balOut.data = bal';
