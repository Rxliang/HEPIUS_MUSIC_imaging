function harmonicIndex = getHarmonicIndices(env1, env2)
%getHarmonicIndices Rewrite version of TCD3D code to get harmonic indices
%used to generate ICP factor 15

L = length(env1) - 1; 
N = 8;

% Initialize
env1pad(1:L*N) = 0;
env2pad(1:L*N) = 0;

% repeat envelope N times
count = 1;
for idx = 1:N 
    env1pad(count:count + L-1) = (env1(1:L));
    env2pad(count:count + L-1) = (env2(1:L));
    count = count + L;
end

dt = 1/L;
fftLen = 10*L;
fftN = 1:fftLen;
fftN2 = fftLen/2;
dfft = 1/(fftLen*dt);
DaznioAsis = fftN.*dfft;

f1 = 0.05/dfft; 
f2 = 0.5/dfft;
f3 = 0.8/dfft; 
f4 = 1.2/dfft;

% create fft's
env1spec = abs(fft(env1pad(:)));
env2spec = abs(fft(env2pad(:)));

ratio2_1 = env2spec(9)/env2spec(1);
ratio1_1 = env1spec(9)/ env1spec(1); 

harmonicIndex.index1 = ratio2_1./ratio1_1; 

ratio2_4 = [env2spec(9)/env2spec(1) env2spec(17)/env2spec(9) env2spec(25)/env2spec(33)] ;
ratio1_4 = [env1spec(9)/env1spec(1) env1spec(17)/env1spec(9) env1spec(25)/env1spec(33)];

harmonicIndex.index2 = (sum((ratio2_4-ratio1_4).^2))^0.5;

ratio2_2 = (env2spec(9)^2 + env2spec(17)^2 + env2spec(25)^2 + env2spec(33)^2)^0.5/env2spec(1);
ratio1_2 = (env1spec(9)^2 + env1spec(17)^2 + env1spec(25)^2 + env1spec(33)^2)^0.5/ env1spec(1); 

harmonicIndex.index3 = ratio2_2./ratio1_2; 

ratio2_3 = ([env2spec(9) env2spec(17) env2spec(25) env2spec(33)] )/env2spec(1);
ratio1_3 = ([env1spec(9) env1spec(17) env1spec(25) env1spec(33)])/ env1spec(1); 

harmonicIndex.index4 = (sum((ratio2_3 - ratio1_3).^2))^0.5; 



end

