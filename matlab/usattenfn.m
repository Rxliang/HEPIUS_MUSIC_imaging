function [atten_dB, comp] = usattenfn(depthAx_mm, freq_MHz, alpha)

if nargin < 3
  alpha = 0.54; % dB/MHz
end

atten_dB = alpha * depthAx_mm/10 * freq_MHz;

atten_field = 10.^(atten_dB/20);

comp = atten_field/atten_field(1);

