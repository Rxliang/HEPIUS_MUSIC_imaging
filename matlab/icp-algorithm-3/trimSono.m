function sonoImg = trimSono(sonoImg, threshold, threshType)

% This function is used to trim a fixed number of rows or a
% of rows from the bottom of the sonogram

% It assumes an existing sono for internal, external and for
% both spectrum types (fft and AR)
% threshType = 1: Fixed (default if not specified)
% threshType = 2: Percentage
% TODO Error checking

% Check to see if this is the first time trimming is done
if ~isfield(sonoImg.int, 'fftOriginal')
    % make backups of originals
    sonoImg.int.fftOriginal = sonoImg.int.fft;
    sonoImg.ext.fftOriginal = sonoImg.ext.fft;
    sonoImg.int.AROriginal = sonoImg.int.AR;
    sonoImg.ext.AROriginal = sonoImg.ext.AR;
    return;
end

if nargin < 3,
    threshType = 'fixed';
end

totalRows = size(sonoImg.int.fftOriginal, 1);
switch threshType
    case 'fixed'
        numRowsToCut = threshold;
    case 'prctg'
        numRowsToCut = round(threshold*totalRows/100);
end

sonoImg.int.fft = sonoImg.int.fftOriginal(numRowsToCut + 1:end,:);
sonoImg.ext.fft = sonoImg.ext.fftOriginal(numRowsToCut + 1:end,:);
sonoImg.int.AR = sonoImg.int.AROriginal(numRowsToCut + 1:end,:);
sonoImg.ext.AR = sonoImg.ext.AROriginal(numRowsToCut + 1:end,:);