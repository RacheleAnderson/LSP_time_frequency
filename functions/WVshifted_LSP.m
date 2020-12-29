function [WVshift, TI, FI,WV0]= WVshifted_LSP(noise_level, a_q, b_q, c_q, c_r, dataN, f0, fs, nfft);

% Compute Wigner-Ville spectrum and shift it to f0 centre frequency 

% The Wigner-Ville spectrum of an LSP is W_E (t,theta) = q(t)*R(theta) 
% where R(theta) is the Fourier transform of the function r
% NB. The Wigner distribution is time-shift and frequency-shift invariant

% q(x) = L + a_q* exp(-(c_q/2) * ((x-b_q)^2))
% r(x) =  exp(-(c_r/8).*x.^2)
% R(theta) = 2* sqrt(2*pi / c_r) * exp(- 2/c_r * theta.^2)

tl=[0:dataN-1]'/fs;
f=[-nfft/2:nfft/2-1]'*fs/nfft;

theta = 2*pi*f*ones(1,length(tl));
tau =  ones(length(theta),1)*tl';

Fr = 2 * sqrt(2*pi/c_r) * exp(-(2/c_r)*theta.^2);
q = noise_level + a_q * exp(-c_q/2 .* (tau-b_q).^2);

WV0 = q .* Fr;

% shift VW to match frequency in the signal 
d = (f0)/fs;
L = floor(d*nfft); % the centre freq is at the L-th row
WVshift  = WV0(floor(nfft/2-L):nfft/2+floor(nfft/2-L)-1,:); 


TI=tl;
FI=[0:nfft/2-1]'/nfft*fs; 

WVshift = WVshift'/2;  


