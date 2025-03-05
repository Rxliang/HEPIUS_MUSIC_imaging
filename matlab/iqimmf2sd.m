% mutiframe iq version

vtData = getenv('VITTAMED_DATA');
matInPath = [vtData '/nicp/iqdata/'];
startTime_s=0;
if 1
  inFilePre = '_2_20180710164348365';
  inFilePre = '_1_20180710220239819';
  inFilePre = '_2_20180710220236785';
  inFilePre = 'site1_subject1_20180711000451318_1';
  
  inFilePre = 'bcm_test_20180713c_20180713143438668_2'
end

inFile = [matInPath inFilePre]

iqSet = readreimfn(inFile);

IQVec = [];
for k = 1:size(iqSet.IQMat,3)
  IQVec(k) = sumall(iqSet.IQMat(:,:,k));
end
    
spec.HPFOrder = 6;
spec.cutHPF_Hz = 100;
spec.msWinLen = 15;
spec.winOffsetFrac = 0.5;
spec.spectralMethod = 'burg'; % fft, yule
spec.binThresh = 1e-1;

plt=1;
figure(1)
sonoSession1=[];
env1=[];
spec.T_s = iqSet.tAx(end)-startTime_s;
[sonoSession1, env1, ...
   sumPowerEnv1] = calcsonofn(IQVec, spec, plt);

return


pwrThres = 0.035;
DopState = 'computeCFIPowerEst';
 
im=[];
clear cdiprocoffline

q=1;
np=20;
cnti=1;
pwr = [];
tAx_s = 0;
frameDelay_s = np*ttnaDop_s;

mn = [];
mx = [];

figure(1)
clf 
while q+np-1 <= size(IQMat,3)
  qStart = q;
  qEnd = qStart+np-1;
  iqm = IQMat(:,:,qStart:qEnd);
  q=qEnd+1;

  sziqm = size(iqm);
  
  iqmThis = reshape(iqm, [sziqm(1:2) 1 np]); 
  
  %  iqmThis = iqmThis(15:20, 7:12, :, :);
  
  im(:,:,cnti) = abs(cdiprocoffline(iqmThis));

  %  pwr(cnti) = sumall(abs(im(:,:,cnti)));
  pwr(cnti) = sumall((im(:,:,cnti)));
  
  mn(cnti) = minall(im(:,:,cnti));
  mx(cnti) = maxall(im(:,:,cnti));  
  
  cnti=cnti+1;
  tAx_s(cnti) = tAx_s(cnti-1)+frameDelay_s;
  
end

tAx_s = tAx_s(1:length(pwr));
figure(2)
clf
plot(tAx_s,pwr)

winLow = median(mn);
winHigh = median(mx);

if 1
figure(1)
clf
colormap gray
for q = 1:size(im,3)
  if 1
      imagesc(im(:,:,q), [winLow winHigh]);
      drawnow
      pause(frameDelay_s);
  end
end
end

if 0
% make B-mode image

sumIQ = IQMat(:,:,1); %sum(IQMat,3);
env = abs(sumIQ);
env_dB=20*log10(env);
env_dB=env_dB-max(max(env_dB));
env_gray=255*(env_dB+60)/60;

figure(3)
imagesc(env_gray);
colormap gray

end

figure(3)
rs=9:12;
cs=5:8;
spec.binThresh = 1e-9;
for r = rs
  for c = cs
    IQVec = IQMat(r,c,:);
    sonoSession1=[];
    env1=[];
    spec.T_s = tAx(end)-startTime_s;
    [sonoSession1, env1, ...
     sumPowerEnv1] = calcsonofn(IQVec, spec, plt);
    pausede
  end
end


rowDecim = 4;
colDecim = 4;

ovlSz = floor([sz(1)/rowDecim, sz(2)/colDecim]);

overlay = zeros(ovlSz);
sumPowerEnv=[];

spec.T_s = tAx(end)-startTime_s;
spec.binThresh = 1e-9;

if  ~exist('sonoSession', 'var') | isempty(sonoSession)
    
    sonoSession=[];
    env = [];
    cntr=1;
    cnts=0;
    for k = 1:rowDecim:sz(1)
        cntc=1;
        for l = 1:colDecim:sz(2)
            rl = k;
            rh = min(k+rowDecim-1,sz(1));
            cl=l;
            ch=min(l+colDecim-1, sz(2));
            
            [rl rh cl ch]
            IQCrop = IQMat(rl:rh, cl:ch, :);        
            IQVec = sum(sum(IQCrop,1),2)/(rowDecim*colDecim);
            plt=1;
            cnts=cnts+1;
            [sonoSession{cnts}, env(cnts,:), ...
             sumPowerEnv(cntr,cntc)] = calcsonofn(IQVec, spec, plt);
            
            %    if cntr==9 & cntc==7
            if cntr==4 & cntc==4
                %pausede
            end
            
            %    pausede
            cntc=cntc+1;
        end
        cntr=cntr+1;
    end
    
end
figure(1)
clf

pltEst=1;
period=[];
FAbs=[];
q=0;
for k = 1:cntr
  for l = 1:cntc
    q=q+1;
    if q > size(env,1)
        break;
    end
    
    inSig = env(q,:);
    fs = 1/(sonoSession{q}.timeAx(2)-sonoSession{q}.timeAx(1)); 
    figure(1)
    clf
    plot(sonoSession{q}.timeAx,inSig);
    figure(2)
    clf
    [period(k,l), FAbsPre, fAxis] = estimateperiod(inSig, fs, fs, ...
                                              pltEst);
    FAbs(q,:) = FAbsPre.';
    
    q
    pausede

  end
end

spec.binThresh=1e-10;
for r = 1:cnts
    imgInt1 =  sonoSession{r}.sonoWav.int;
    sonoPlot=[];
    sonoPlot(1).isLog = 1; 
    sonoPlot(1).isGray = 1; 
    sonoPlot(1).title = 'Low quality image'; 
    sonoPlot(1).img = imgInt1;
    sonoPlot(1).env = [];
    imBin = sonoSession{r}.binarizeSono(imgInt1, spec.binThresh); % %.0005);
    y3 = sonoSession{r}.getEnvelopeFromBinaryImg(imBin,1);
    %y3 = TVL1denoise(y3, 1.5, 100);
    %y4a = medfilt1(y4,3);
    %sonoPlot(1).img = imBin; %imgExt1;
    sonoPlot(1).env = y3;


    figure(1);
    sonoSession{r}.plotSono(sonoPlot);
    pausede
end





  return

    


IQMat = squeeze(iqr(r,c,:)+i*iqi(r,c,:));
szIQMat = size(IQMat);

IQVec = sum(sum(IQMat,1),2)/prod(szIQMat(1:2));


spec.T_s = framePeriod_s*numFrames-startTime_s;





return


sonoSession = sono;

if 1
  if 0
    % good quality settings
    sonoSession.msWinLen = 50;
    sonoSession.winOffsetFrac = 0.05;
  end
  
  if 2
    % economy settings
    sonoSession.msWinLen = 100;
    sonoSession.winOffsetFrac = 0.5;
  end

  sonoSession.rawIQ(:,1) =  real(IQVec);
  sonoSession.rawIQ(:,2) =  imag(IQVec);
  sonoSession.T_s = framePeriod_s*numFrames-startTime_s;
  sig=1;
else
  a=1;
  b=1;
  sets{a}.adcFiles{b} = 'ADC_024.dat';
  dataPath = '~/Documents/icp-algorithm/data/aarau/data/2015_03_15__15_52/Data/';
  sonoSession.loadAdcFile(sets{a}.adcFiles{b}, dataPath);
  sonoSession.T_s = 10;
  sonoSession.msWinLen = 50;
  sonoSession.winOffsetFrac = 0.05;
  sig=1;
end


sonoSession.HPfilterOrder = 6; %2 HPF order
sonoSession.cutHPF_Hz = 100; % 100 HPF cutoff
        
sonoImg.int.fft = [];
%sonoSession.gensonofwdrevfn(1, 'burg');
%sonoSession.gensonofwdrevfn(1, 'fft');

tic
%sonoSession.gensonofwdrevfn(sig, 'yule');
sonoSession.gensonofwdrevfn(1, 'burg')
%sonoSession.gensonofwdrevfn(1, 'fft');
toc

if sig==1
  sonoImg.int.AR = sonoSession.sonoWav.int;
else
  sonoImg.int.AR = sonoSession.sonoWav.ext;
end

%imagesc(sonoSession.timeAx, sonoSession.fAxSigned, ...
%        flipud(log(sonoImg.int.AR)), [45 62])
%colormap(gray)


sonoSum = sum(sonoSession.sonoWav.int,1);

imgInt1 =  sonoSession.sonoWav.int;

sonoPlot(1).isLog = 1; 
sonoPlot(1).isGray = 1; 
sonoPlot(1).title = 'Low quality image'; 
sonoPlot(1).img = imgInt1;
sonoPlot(1).env = [];
imBin = sonoSession.binarizeSono(imgInt1, 1e-9); % %.0005);
y3 = sonoSession.getEnvelopeFromBinaryImg(imBin,1);
%y3 = TVL1denoise(y3, 1.5, 100);
%y4a = medfilt1(y4,3);
%sonoPlot(1).img = imBin; %imgExt1;
sonoPlot(1).env = y3;


figure(1);
sonoSession.plotSono(sonoPlot);
