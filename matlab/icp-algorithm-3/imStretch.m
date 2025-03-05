function imgOut = imStretch(img)
%inStretch Recasts sono img to be between 0 and 1
%   Detailed explanation goes here

finf = isinf(img);

minImg = min(min(img(~finf)));

imgOut = img - minImg;
maxImg = max(max(imgOut(~finf)));
imgOut = imgOut./maxImg;

imgOut(finf) = 1;

end

