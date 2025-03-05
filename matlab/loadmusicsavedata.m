vtData = getenv('MUSIC_DATA');
dataPath = fullfile(vtData);
%outFilePre = fullfile(dataPath, mfile);

% these are limited frame snapshots, not cine

fn = '20231214_MUSIC_Plane_Sag_WeakFlowVLow_Tubing.mat';
fn = '20231214_MUSIC_Plane_Axial_WeakVFlowLow.mat';
%fn = '20231214_MUSIC_FocusedDoppler_WeakFlowVLow.mat';
inMat = fullfile(dataPath, 'pig', ...
                 fn);

%if ~exist('d', 'var') | isempty(d)
%  d = load(inMat);
%end

if ~exist('Receive', 'var') | isempty(Receive)
    %  d = load(inMat);
  load(inMat);
end

% for Plane_Axial

ev = [];
ev.tx = [Event.tx];

% 51:61 are pulse dop with diff Doppler angles
% 56 os zero degree
ind = find(ev.tx == 56);

indDiff = find(diff(ind)~=11);

indEndEns = ind(indDiff);
indEndEns = [1 indEndEns ind(end)];

numEns = length(indEndEns);

indEns = [];
evEns = [];
rcvEns = [];
for q = 2:numEns
  indF = find( (ind <= indEndEns(q)) & (ind > indEndEns(q-1)) )
  indEns{q-1} = ind(indF);
  evEns{q-1} = Event(indEns{q-1});
  rcvEns{q-1} = [Event(indEns{q-1}).rcv];
end

RF = [];

elSel = 1:64;
ttna_us = [];
seqControlDop = 6;
ttnaFix_us = ceil(1/(dop.PRF*dop.numAngs)*1E6);
%ttnaFix_us = 84; % needs to propogate via mat file in future

for q = 1:numEns-1
    rcvThis = rcvEns{q};
    evThis = evEns{q};
  for r = 1:length(rcvThis)
    seqc = evThis(r).seqControl;
    ttna_us(q,r) = ttnaFix_us;
    thisReceive = Receive(rcvThis(r));
    st = thisReceive.startSample;
    en = thisReceive.endSample;
    sd = thisReceive.startDepth;
    ed = thisReceive.endDepth;    
    dmFreq_MHz = thisReceive.demodFrequency;
    lambda_mm = 1e-3*P.c/dmFreq_MHz;
    dz_mm = lambda_mm/thisReceive.samplesPerWave;
    bufNo = thisReceive.bufnum;
    frameNo = thisReceive.framenum;
    RF(:,:, r, frameNo) = RcvData{bufNo}(st:en, elSel, frameNo);      
  end
end


lambda_mm = 1e-3*P.c/dmFreq_MHz;
rfIm = log(squeeze(abs(RF(:,32,:,4))));
ensInt_us = dop.numAngs/dop.PRF;

dzRF_mm = (ed-sd)*lambda_mm/(en-st+1);
dtRF_us = ttnaFix_us*dop.numAngs;

axRF = makeimaxisfn(rfIm, [dtRF_us dzRF_mm]);
axRF.vAxB = axRF.vAxB - axRF.vAxB(1);
axRF.vAxC = axRF.vAxC - axRF.vAxC(1);

for f = 1:4
figure(1)
imagesc(axRF.uAxB, axRF.vAxB, log(squeeze(abs(RF(:,32,:,f)))))
xlabel('\mus ');
ylabel('depth (mm)');
pausede
end



return



% ImgData1: {170×223×1×2 double}
                                  %wid      %ht in lambda
% Pdata1  Size: [170 223 1], delta [0.4566 0 1]
lambda_mm = 1e-3*P.c/dmFreq_MHz;
ax = makeimaxisfn(ImgData{1}(:,:,1,1), lambda_mm*PData(1).PDelta([1 3]));
ax.vAxB = ax.vAxB-ax.vAxB(1);
ax.vAxC = ax.vAxC-ax.vAxC(1);

for f = 1:10
figure(1)
clf
%imagesc(log(d.ImgData{1}(:,:,1)))%imagesc(log(d.ImgData{1}(:,:,1)))
%imagesc(ax.uAxB, ax.vAxB, log(ImgData{1}(:,:,1,f)))
imagesc(log(ImgData{2}(:,:,1,f)))
colormap gray
xlabel('mm');
ylabel('mm');
drawnow
pausede
end


axP = makeimaxisfn(ImgData{2}(:,:,1,1), lambda_mm*PData(2).PDelta([1 3]));
axP.vAxB = axP.vAxB-axP.vAxB(1);
axP.vAxC = axP.vAxC-axP.vAxC(1);

figure(2)
clf
%imagesc(log(d.ImgData{1}(:,:,1)))%imagesc(log(d.ImgData{1}(:,:,1)))
imagesc(axP.uAxB, axP.vAxB, (ImgData{2}(:,:,1,5)))
colormap gray
xlabel('mm');
ylabel('mm');
