mfile = mfilename;
addpath('icp-algorithm-3')

spec = [];
spec.HPFOrder=6;
spec.cutHPF_Hz=50;
spec.msWinLen = 50;
spec.winOffsetFrac = 0.05;
spec.spectralMethod = 'burg';   spec.binThresh=8e-14;
%  spec.spectralMethod = 'fft';   spec.binThresh=9e-7;

  
spec.T_s(1) = dop.numAngs*iqSet.ttna_s;
plt=1;
iqVecInt = sum(iqSet.IQMat, 2);
iqVecInt = squeeze(sum(iqVecInt, 1)).';
  
iqVec = [];
iqVec{1}  = iqVecInt;  

[sonoSession, env, sumPowerEnv] = calcsonofn(iqVec, spec, plt)
  
  
%  sono = gensonofwdrevfn(I, Q, fs, params, plt)
 


