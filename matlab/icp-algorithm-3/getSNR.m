function snr = getSNR(env, img, snrMethod)

% computes snr of envelope and or img according to strategy set in method

% power spectrum envelope method
envNorm = env./prctile(env,95); 

% filter to smooth
[B, A] = butter(4, .1, 'low');
envFilt = filtfilt(B,A, envNorm);

% difference between smoothed and unsmoothed
dEnv = envNorm - envFilt;

% Welch power spectrum of diff
nfft = 64;
[pxx, f] = pwelch(detrend(dEnv), hamming(nfft), round(nfft*.8), nfft);

% energy in high freq side of spectrum
snr.energy = sqrt(sum(pxx(round(nfft/4)+1:end)));