mfile = mfilename;
%definebnsenv;
vtHome = getenv('MUSIC_HOME');
vtData = getenv('MUSIC_DATA');
inPath = fullfile(vtData, 'nicp', 'iqdata');
inFilePre = '_20190223010517346_1'; % blank
inFilePre = '_20190223010517286_2'; % blank
inFilePre = '_20190223132424786_1'
inFilePre = '_20190223132429800_2' % short

%inFilePre = '_20190223132409618_2'; % poor

inFilePre = '_20190223131955240_1'; 

startTime_s = 0;

[outFile, iqSet] = readiqv2fn(inFilePre, startTime_s, inPath, 2 , 1);


  spec.HPFOrder=6;
  spec.cutHPF_Hz=50;
  spec.msWinLen = 50;
  spec.winOffsetFrac = 0.05;
  spec.spectralMethod = 'burg';   spec.binThresh=8e-14;
%  spec.spectralMethod = 'fft';   spec.binThresh=9e-7;

  
%spec.T_s(1) = size(iqSet{s,1}.IQMat,3)*iqSet{s,1}.ttna_s;
% spec.T_s(2) = size(iqSet{s,2}.IQMat,3)*iqSet{s,2}.ttna_s;  
  spec.T_s(1) = size(iqSet.IQMat,3)*iqSet.ttna_s;
  spec.T_s(2) = size(iqSet.IQMat,3)*iqSet.ttna_s;  
  
  plt=1;
  %iqVecInt = sum(iqSet{s,1}.IQMat, 2);
  iqVecInt = sum(iqSet.IQMat, 2);
  iqVecInt = squeeze(sum(iqVecInt, 1)).';
  %  iqVecExt = sum(iqSet{s,2}.IQMat, 2);
  %iqVecExt = squeeze(sum(iqVecExt, 1)).';  
  
  iqVec = [];
  iqVec{1}  = iqVecInt;
  iqVec{2}  = iqVecInt;
  
  [sonoSession, env, sumPowerEnv] = calcsonofn(iqVec, spec, plt)
  
  