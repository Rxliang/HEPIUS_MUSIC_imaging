function paddedSignal = zeropad(inSignal)

% Appends 0 along dimension of vector so total length is power of 2

[m, n] = size(inSignal);

if m > 1
    if mod(log2(m),1) == 0
        paddedSignal = inSignal;
    else
        L = floor(log2(m))+1;
        paddedSignal = [inSignal; zeros(2^L-m,1)];
    end
elseif n > 1
    if mod(log2(n),1) == 0
        paddedSignal = inSignal;
    else
        L = floor(log2(n))+1;
        paddedSignal = [inSignal zeros(1,2^L-n)];
    end
else
    % data is a scalar and should throw and error
    error('Data is a scalar');
end