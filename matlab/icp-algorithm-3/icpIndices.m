function icpIndex = icpIndices(env1, env2)

% Function to take an envelope as input and generate indices for
% calculating ICP balance as per TCD3D factoirs 2, 5 , 12 and 15

% find enevlope above 70% of maximum
env1_70 = find(env1 > 0.7*max(env1));
env2_70 = find(env2 > 0.7*max(env2));
% normalize by mean of these values
env1_70norm = env1/mean(env1(env1_70));
env2_70norm = env2/mean(env2(env2_70));

% fit a 1st order polynomial between env1 and env2
[p, s] = polyfit(env1, env2, 1); 
% find the value of the polynomial at the max of env1
polyNormFactor = polyval(p, max(env1));
        
% Normalize envelopes by max value
env1_maxNorm = env1/max(env1);
env2_maxNorm = env2/max(env2);

% Normalize env 2 by polyNormFactor
env2_polyNorm = env2/polyNormFactor(1);

% Create approximate linear fit to normalized envelopes
N = 20;
envMin = min(env1_maxNorm);
envMax = max(env1_maxNorm);
dx = (envMax - envMin)/N;
count = 0;
Env1_fit = [];
Env2_fit = [];
for idx = 1:N-1
    a1 = envMin + dx*(idx-1);
    a2 = envMin + dx*idx;
    f = find(env1_maxNorm >= a1 & env1_maxNorm < a2);
    
    if ~isempty(f)
        count = count+1;
        Env1_fit(count) = mean(env1_maxNorm(f));
        Env2_fit(count) = mean(env2_polyNorm(f));
    end
end

% generate indices
icpIndex.index5N = mean(abs(env1_maxNorm - env2_polyNorm)); %% 5
icpIndex.index5 = mean(abs(env1_maxNorm - env2_maxNorm)); %% 5

icpIndex.index2 = std(env1_maxNorm - env2_polyNorm); % 2

icpIndex.index12 = mean(abs(Env1_fit - Env2_fit)); % 12

hIndex = getHarmonicIndices(env1, env2);
icpIndex.index15 = hIndex.index4;

