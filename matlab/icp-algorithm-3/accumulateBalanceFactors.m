function balSet = accumulateBalanceFactors(bal)

% Takes a balance factor structure and accumulates 
[r, c] = size(bal);

% get max length
Lmax = 0;

for a = 1:r
    PIenv = [];
    for b = 1:c
        if bal{a,b}.snr > 0
            PIenv = [PIenv bal{a,b}.PIenv.all];
        end
        L = length(PIenv);
        if L > Lmax
            Lmax = L;
        end
    end
end

balSet.PIenv = ones(r, L)*NaN;
balSet.mom1 = ones(r, L)*NaN;
balSet.mom2 = ones(r, L)*NaN;
balSet.hrEnv = ones(r, L)*NaN;
balSet.hrMom1 = ones(r, L)*NaN;
 
for a = 1:r
    PIenv = [];
    mom1 = [];
    mom2 = [];
    hrEnv = [];
    hrMom1 = [];
    for b = 1:c
        if bal{a,b}.snr > 0
            PIenv = [PIenv bal{a,b}.PIenv.all];
            mom1 = [mom1 bal{a,b}.PImom1.all];
            mom2 = [mom2 bal{a,b}.PImom2.all];
            hrEnv = [hrEnv bal{a,b}.hrEnv.all];
            hrMom1 = [hrMom1 bal{a,b}.hrMom1.all];
        end
        
    end
    balSet.PIenv(a, 1:length(PIenv)) = PIenv;
    balSet.mom1(a, 1:length(PIenv)) = mom1;
    balSet.mom2(a, 1:length(PIenv)) = mom2;
    balSet.hrEnv(a, 1:length(PIenv)) = hrEnv;
    balSet.hrMom1(a, 1:length(PIenv)) = hrMom1;
end

% augment missing as nan
f = find(balSet.PIenv == 0);
balSet.PIenv(f) = NaN;
f = find(balSet.mom1 == 0);
balSet.mom1(f) = NaN;
f = find(balSet.mom2 == 0);
balSet.mom2(f) = NaN;
f = find(balSet.hrEnv == 0);
balSet.hrEnv(f) = NaN;
f = find(balSet.hrMom1 == 0);
balSet.hrMom1(f) = NaN;

