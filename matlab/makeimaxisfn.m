function [ax] = makeimaxisfn(im, pixWid)

sz = size(im);

uLimB = sz(2)/2 * pixWid(1);
vLimB = sz(1)/2 * pixWid(2);

ax.uAxB = -uLimB:pixWid(1):uLimB;
ax.vAxB = -vLimB:pixWid(2):vLimB;

ax.uAxC = axbound2centre(ax.uAxB);
ax.vAxC = axbound2centre(ax.vAxB);
