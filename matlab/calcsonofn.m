function [sonoSession, env, sumPowerEnv] = calcsonofn(IQVec, spec, plt)

sonoSession = sono;

sonoSession.rawIQ{1}(:,1) =  real(IQVec{1});
sonoSession.rawIQ{1}(:,2) =  imag(IQVec{1});

for sig = 1:1
  %sig=2;
  sonoSession.HPfilterOrder = spec.HPFOrder;
  sonoSession.cutHPF_Hz = spec.cutHPF_Hz;
  sonoSession.msWinLen = spec.msWinLen;
  sonoSession.winOffsetFrac = spec.winOffsetFrac;
  sonoSession.posPGramFlag = spec.posPGramFlag; % show +ve freqs only
  sonoImg.int.fft = [];
  sonoSession.T_s = spec.T_s;

  tic
    sonoSession.gensonofwdrevfn(sig, spec.spectralMethod);
    %sonoSession.gensonofwdrevfn(sig, 'yule');
    %sonoSession.gensonofwdrevfn(1, spec.spectralMethod);
    %sonoSession.gensonofwdrevfn(1, 'fft');
  toc

  if sig==1
    sonoImg.int.AR = sonoSession.sonoWav.int;
    sonoSession.fAx = sonoSession.fAx.int;
    sonoSession.fAxSigned = sonoSession.fAxSigned.int;  
    sonoSession.timeAx = sonoSession.timeAx.int;  
    img =  sonoSession.sonoWav.int;
  else
    sonoImg.ext.AR = sonoSession.sonoWav.ext;
    sonoSession.fAx = sonoSession.fAx.ext;
    sonoSession.fAxSigned = sonoSession.fAxSigned.ext;    
    sonoSession.timeAx = sonoSession.timeAx.ext;
    keyboard
    img = sonoSession.sonoWav.ext;
  end

  imBin = sonoSession.binarizeSono(img, spec.binThresh); 
  imMasked = img.*imBin;
  ind = find(sonoSession.fAx < sonoSession.cutHPF_Hz);
  imMasked(ind,:)=0;

  env = sonoSession.getEnvelopeFromBinaryImg(imBin,1);

  for i = 1:size(imMasked,2)
    r = round(env(i));
    if r < 0 | r > size(imMasked,2)
      continue
    end
  
    imMasked(round(env(i)):end, i) = 0;
  
  end
  
  sumPowerEnv = sum(imMasked(:));
  
  if plt
    sonoPlot(sig).isLog = 1; 
    sonoPlot(sig).isGray = 1; 
    sonoPlot(sig).title = 'Low quality image'; 
    sonoPlot(sig).img = img;
    sonoPlot(sig).env = [];
 
    %y3 = TVL1denoise(y3, 1.5, 100);
    %y4a = medfilt1(y4,3);
    %sonoPlot(1).img = imBin; %imgExt1;
    sonoPlot(sig).env = env;
    
    figure(sig);
    sonoSession.plotSono(sonoPlot);
  end

end
