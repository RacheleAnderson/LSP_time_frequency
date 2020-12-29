function S  =  scalogram(y,fs,a,wav_name)

% CREDITS: Implementation of the scalogram function is done according to the
% following book:
% Stark, H.G., "Wavelets and Signal Processing. An Application-Based Introduction.",
% Springer-Verlag Berlin Heidelberg, 2005. ISBN 9783540234333

% INPUT: 
%   y:           signal
%   fs:          sampling frequency
%   scale_vec:   vector of scale factors
%   wav_name:    wavelet name

% OUTPUT:
%   scalogram:   scalogram

N = max(size(y));
if size(y,1)>size(y,2)
    y = y'; 
end

%Sampling circular frequency
oms = 2*pi*fs;

%Generation of symmetric time vectors
if floor(N/2) == (N/2)
    t = (-floor(N/2):(floor(N/2)-1))/fs; 
else 
    t = (-floor(N/2):(floor(N/2)))/fs; 
end

rows = length(a); % number of scale factors

% Fourier transforms
yhat = fft(y);

% Matrix-Initialization
matrix = zeros(rows,N);

% Loop for increasing scale factors
for i = 1:rows
    
      psi_scale = conj(feval(wav_name,-t./a(i)));
      psi_scale_hat = fft(psi_scale); % Fourier transform of wavelet transform
    
      % Time translation such that minimal time = 0
      trans = exp((-j*t(1)*(0:(N-1))*oms/N));
      conv_hat = ((yhat.*psi_scale_hat).*trans)/sqrt(a(i));

      matrix(i,:) = ifft(conv_hat);
end

S  =  abs(matrix).^2;