function y = morlet_wav(x)
% Generates the Morlet wavelet

f = pi*sqrt(2/log(2)); % The parameter that allows trade between time and frequency resolutions
y = pi^(-.25)*(exp(-(j*f*x))-exp(-(f*f)/2)).*exp(-(x.*x)/2);