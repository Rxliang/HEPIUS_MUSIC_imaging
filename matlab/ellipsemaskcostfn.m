function [cost, ellIm] = ellipsemaskcostfn(prm, im, X, Z, tuneParms)

r0 = prm(1:2);
a = prm(3);
b = prm(4);
artVal = prm(5);

[ellImPre, discInd] = drawellipsefn(r0, ...
                                a, ...
                                b, X, Z);
ellIm = artVal*ellImPre;

%cost = sum((ellIm(discInd) - im(discInd)).^2);
cost = sum((ellIm(:) - im(:)).^2);

