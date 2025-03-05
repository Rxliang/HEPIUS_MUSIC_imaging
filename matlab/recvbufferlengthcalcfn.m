function [samplesPerEvent, maxDist_wvl] = recvbufferlengthcalcfn(FOV)

lambda_m = FOV.c_mps/FOV.freq_MHz/1e6;
samplesPerWvlThis = ceil(FOV.baseSamplesPerWvl * FOV.freq_MHz/FOV.baseFreq_MHz);

maxDistanceForRcvEndDepth_mm = sqrt(FOV.maxDepth_mm.^2 + (FOV.maxHalfWid_mm*2).^2);

maxDistanceForRcvBufferCalc_mm = ...
    maxDistanceForRcvEndDepth_mm - FOV.startDepth_mm;

lambdaBase_mm = FOV.c_mps/FOV.baseFreq_MHz*1e-3;

maxDist_wvl_Pre = maxDistanceForRcvBufferCalc_mm/lambdaBase_mm;
numberOfWvlRoundTrip = maxDist_wvl_Pre*2;

samplesPerEventPre = numberOfWvlRoundTrip* ...
                                     samplesPerWvlThis;

% receive buffer always uses multiples of 128
samplesPerEvent = findnextmultiplefn(samplesPerEventPre, 128);
samplesPerEventAdded = samplesPerEvent-samplesPerEventPre;

% for Receive.endDepth, just use caclulated value

maxDist_wvl = maxDistanceForRcvEndDepth_mm/lambdaBase_mm + samplesPerEventAdded/samplesPerWvlThis/2;
