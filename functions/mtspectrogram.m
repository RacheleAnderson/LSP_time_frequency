function [S,TI,FI]=mtspectrogram(X,WIN,Fs,NFFT,NSTEP,wei);

% MTSPECTROGRAM MTSpectrogram 
% [S,TI,FI]=mtspectrogram(X,WIN,Fs,NFFT,NSTEP,wei) calculates the spectrogram using the 
% short-time Fourier transform of the (NN X 1) data vector X and a
% window WIN. A length of a normalized Hanning window can be specified
% as WIN. Any other window function can be also set as the parameter
% WIN, or alternatively a set of multitapers as a matrix with a chosen weighting vector, wei.  
%
% INPUT:
%    X:    Data sequence, must be of shorter length than NFFT. 
%    WIN:  A data window of length L less than NN. If the length is 
%          specified as the parameter WIN=L, a Hanning window of length L is used. A matrix
%          of size (L X K) specifiy a set of multitapers. 
%    Fs:   Sample frequency, default Fs=1. 
%    NFFT: The number of FFT-samples, default NFFT=1024.
%    NSTEP:Sample step to the next frame, default=1. 
%    wei:  Weights of the different multitaper spectrograms, size (K X 1), default wei(k)=1/K, k=1 ... K.
%
% OUTPUT:
%    S:    Output spectrogram, matrix of sixe (NN X NFFT).
%    TI:   Time vector (NN X 1).
%    FI:   Frequency vector (NFFT X 1).

if nargin<2
    'Error: No data input and window'
end
[L,K]=size(WIN);
if L<K
  [L,K]=size(WIN');
  WIN=WIN';
end
if L==1 & K==1
    L=WIN;
    WIN=hanning(WIN);
    WIN=WIN./sqrt(WIN'*WIN);
end
if nargin<6
    wei=ones(K,1)/K;
end

if nargin<5
    NSTEP=1;
end

if nargin<4
    NFFT=1024;
end
if nargin<3
    Fs=1;
end

[NN,M]=size(X);
if M>NN
    X=transpose(X);
    NN=M;
end
wei=wei(:);

x=[zeros(fix(L/2),1);X;zeros(ceil(L/2),1)];
S=[];
TI=[];

for i=1:NSTEP:length(x)-L
  testdata=x(i+1:i+L);
  F1=abs(fft(WIN.*(testdata*ones(1,K)),NFFT)).^2;
  Fsave=F1*wei;
  S=[S;Fsave(1:NFFT)'];
  TI=[TI;i];
end

TI=TI/Fs;
FI=[0:NFFT/2-1]'/(NFFT)*Fs;

% Correcting for sampling frequency
S=S(:,1:NFFT/2)/Fs;

% figure
%%% S=fftshift(S(1:NFFT,:),1);
% c=[min(min(S)) max(max(S))];
% pcolor(TI,FI,S')  
% shading interp
% caxis(c)
% ylabel('Frequency (Hz)')
% xlabel('Time (s)')
% title('Spectrogram')