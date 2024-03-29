
% This script allows to reproduce the simulation study presented in the paper:
% Anderson, R., Sandsten, M. "Time-frequency feature extraction for classification of episodic memory"
% EURASIP J. Adv. Signal Process. 2020, 19 (2020). https://doi.org/10.1186/s13634-020-00681-8

% -------------------------------------------------------------------------

% The code in MSE_methods_comparison.m does the following:
% 1 - Simulate realizations 
% 2 - Compute exact Wigner-Ville spectrum (WVS)
% 3 - Evaluate mean square error (MSE) with respect to exact WVS using the following state-of-the-art estimators:
%	3.i     - Hanning window spectrogram (HANN)
%	3.ii    - Welch method with 50% window overlap (WOSA)
%   3.iii   - classical Wigner-Wille spectrum estimate (WV)   
%	3.iv    - optimal LSP kernel with true parameters (LSP)
%	3.v     - optimal LSP kernel with estimated	parameters (LSP-HATS)
%   3.vi    - Continuous wavelet transform with Morlet wavelet (CWT) 

% Each method is optimized to evaluate it at its best performance in terms of MSE:
% - for HANN we optimize the window length
% - for WOSA we optimize the number of windows
% - for WV and CWV we adjust the estimated magnitude with an optimized scaling
% parameter

% Matlab toolboxes used: Symbolic Math Toolbox, Wavelet Toolbox 
% -------------------------------------------------------------------------
% Copyright Anderson, R., Sandsten, M.

clear all; close all
addpath('functions')

%% 1 - Simulate LSP realizations 

num_real = 100; % number of realizations % change to a smaller number (e.g. 10) for quicker results

rng(50) % set seed for reproducibility 

% Model parameters:
a_q = 500; b_q = 0.2; c_q = 800; c_r = 15000; noise = 120; % model parameters
dataN = 256; % samples in each realization
T0 = 0; % initial time
Tf = 0.5; % final time
delta_t = abs(Tf-T0)/(dataN-1); % sampling interval
time_vec = T0 + [0:dataN-1]'* delta_t; % vector of times

fs= 512;
f0 = 25;  % centre frequency f0 (Hz)

% simulate realizations:
[X,X_freq,C,C_freq,R,R_freq,Q] = lsp_f0_sim(num_real,f0,a_q,b_q,c_q,c_r,noise,dataN,time_vec); 

%% 2 - Compute exact Wigner-Ville spectrum (WVS)

nfft = 1024;
[WVSshift,TI,FI,W0] = WVshifted_LSP(noise,a_q,b_q,c_q,c_r,dataN,f0,fs,nfft);
WVS = repmat(WVSshift,1,1,num_real);

%% 3 - Evaluate mean square error (MSE) with respect to exact WVS 

% For the MSE comparison we focus on the area of interest of the time-frequency plane,
% corresponding to 0-100 Hz, [0 - 0.5] seconds
cutSpec = 201; % 100 Hz
WVS = WVS(:,1:cutSpec,:);

%% 3.i - Hanning window spectrogram (HANN)
% For HANN: optimization of the window length

win_length = [dataN/16,dataN/8,dataN/4,dataN/2,dataN]; % window lengths considered
mMSE_HANN_matrix = zeros(length(win_length), num_real); % to store meanMSE computed from several realizations 
                                                        % for every window length
                                                        % the  i:th row contains all the MSE for the num_real realizations
S_HANN = zeros(dataN,nfft/2,num_real,length(win_length)); % to store the spectrograms:
                                                          % S_HANN(:,:,i,j) will be the spectrogram for the i:th realization,
                                                          % computed with window length win_length(j)

for i = 1:length(win_length) 

    win = win_length(i); % window length

    for j = 1:num_real 
        y = X_freq(:,j);
        S = mtspectrogram(y,win,fs,nfft);
        S_HANN(:,:,j,i) = S;        
    end

    mMSE_HANN_matrix(i,:) = compute_MSE(S_HANN(:,1:cutSpec,:,i),WVS); 
       
    clear win SE MSE_f MSE  

end

mMSE_HANN_vec = mean(mMSE_HANN_matrix,2); % mean on the realizations = approximated expected value of the MSE for every window length

[mMSE_HANN_opt,I_HANN_opt] = min(mMSE_HANN_vec); % I_spec_opt is the index corresponding to the MSE optimal parameter for HANN
win_HANN_opt = win_length(I_HANN_opt);
std_MSE_HANN = std(mMSE_HANN_matrix,0,2);
std_MSE_HANN_opt = std_MSE_HANN(I_HANN_opt);

%% 3.ii - Welch method with 50% window overlap (WOSA)

win = dataN; % total length of the num_win windows
num_win = [2,4,8,12,16]; % number of windows considered for optimization
mMSE_WOSA_matrix = zeros(length(num_win),num_real); 
S_WOSA = zeros(dataN,nfft/2,num_real,length(num_win));

for i=1:length(num_win)

     K = num_win(i);
    [WOSA_win, WOSA_wei] = welch_wind(win,K); % Welch windows and weights 1/K

     for j=1:num_real 
        y = X_freq(:,j);
        Sw = mtspectrogram(y,WOSA_win,fs,nfft,1,WOSA_wei);
        S_WOSA(:,:,j,i) = Sw;
        clear y Sw
     end
     
     mMSE_WOSA_matrix(i,:) = compute_MSE(S_WOSA(:,1:cutSpec,:,i),WVS);

     clear K SE MSE_f MSE
end

mMSE_WOSA_vec = mean(mMSE_WOSA_matrix,2);
[mMSE_WOSA_opt,I_WOSA_opt] = min(mMSE_WOSA_vec);
win_WOSA_opt = num_win(I_WOSA_opt);
std_MSE_WOSA = std(mMSE_WOSA_matrix,0,2);
std_MSE_WOSA_opt = std_MSE_WOSA(I_WOSA_opt);

%% 3.iii - Classical Wigner-Wille spectrum estimate (WV)   

WV = zeros(dataN,nfft/2,num_real);

for j = 1:num_real 
    y = X_freq(:,j);
    z = hilbert(y);
    W = quadtf(z,'wigner',[],fs,nfft);    
    WV(:,:,j) = W;

    clear y z W
end

% Find optimal spectral magnitude scaling parameter alpha 
WV = WV(:,1:cutSpec,:);
num = sum(sum(mean(WV.*WVS,3))); 
den = sum(sum(mean(WV.^2,3))); 
alpha_WV = num./den;

MSE_vec_WV = compute_MSE(alpha_WV*WV,WVS);
MSE_vec_WV = squeeze(MSE_vec_WV);
mMSE_WV = mean(MSE_vec_WV); 

% std_MSE_WV = std(MSE_vec_WV);

%%	3.iv and 3.v - LSP optimal kernel with true parameters (LSP) and LSP optimal kernel with estimated	parameters (LSP-HATS)
 
% Settings for HATS inference
T0 = time_vec(1); Tf = time_vec(end);
% Lower bounds:
noise_LB = 0; a_q_LB = 100; b_q_LB = T0; c_q_LB = 1; c_r_LB = 1; 
% Upper bounds:
noise_UB = 150; a_q_UB = 2000;  b_q_UB = Tf; c_q_UB = 10000; c_r_UB = 70000; 
% Starting points for optimization
a_q_0 = (a_q_UB-a_q_LB)/2;
b_q_0 = (b_q_UB-b_q_LB)/2;
noise_0 = (noise_LB+noise_UB)/2;
c_q_0 = (c_q_LB+c_q_UB)/2;
c_r_0 = (c_r_LB+c_r_UB)/2;

LB = [noise_LB, a_q_LB, b_q_LB, c_q_LB, c_r_LB];
UB = [noise_UB, a_q_UB, b_q_UB, c_q_UB, c_r_UB];
THETA_0 = [noise_0, a_q_0, b_q_0, c_q_0, c_r_0];

% NB: Use X, without frequency, because HATS doesn't assume centre frequency 
[noise_HATS,a_q_HATS,b_q_HATS,c_q_HATS,c_r_HATS] = HATS(X,time_vec,LB,UB,THETA_0);

% compute LSP windows and weights
[uopt_HATS,sopt_HATS] = optimal_kernel_LSP(dataN,a_q_HATS,c_q_HATS,c_r_HATS,noise_HATS,fs);
[uopt_true,sopt_true] = optimal_kernel_LSP(dataN,a_q, c_q,c_r,noise,fs);

% compute LSP spectrograms
S_LSP_HATS = zeros(dataN,nfft/2,num_real);
S_LSP_true_par = zeros(dataN,nfft/2,num_real);

for j = 1:num_real 
    
    y = X_freq(:,j); % current realization
    S_LSP_HATS_curr = mtspectrogram(y,uopt_HATS,fs,nfft,1,sopt_HATS); 
    S_LSP_true_curr = mtspectrogram(y,uopt_true,fs,nfft,1,sopt_true); 
    
    S_LSP_HATS(:,:,j) = S_LSP_HATS_curr;
    S_LSP_true_par(:,:,j)= S_LSP_true_curr;  
    
    clear S_LSP_HATS_curr S_LSP_true_curr
    
end

% Compute and save mean square error (LSP true par)
MSE_LSP_true_par = compute_MSE(S_LSP_true_par(:,1:cutSpec,:),WVS);
mMSE_LSP_true_par = mean(MSE_LSP_true_par); % mean on the realizations
std_MSE_LSP_true_par = std(MSE_LSP_true_par); % std on the realizations

% Compute and save mean square error (LSP HATS)
MSE_LSP_HATS = compute_MSE(S_LSP_HATS(:,1:cutSpec,:),WVS);
mMSE_LSP_HATS = mean(MSE_LSP_HATS); % mean on the realizations
std_MSE_LSP_HATS = std(MSE_LSP_HATS); % std on the realizations

    
%%   3.vi - Continuous wavelet transform with Morlet wavelet (CWT) 

f_vec = FI(2:cutSpec+1); % from this vector of frequencies we derive the scales for the wavelets
a_vec = centfrq('morl')./f_vec;

CWT = zeros(dataN,cutSpec,num_real);

for j=1:num_real 

    y= X_freq(:,j);   
    Sc = scalogram(y,fs,a_vec,'morlet_wav'); % S is a matrix of size is length(a_vec) x length(y)     
    CWT(:,:,j)= Sc';

    clear y Sc
    
end

% Find optimal alpha
num_CWT = sum(sum(mean(WVS.*CWT,3))); % .* element-wise multiplication 
den_CWT = sum(sum(mean(CWT.^2,3))); 
alpha_CWT = num_CWT./den_CWT;

% Compute and save mean square error 
MSE_vec_CWT = compute_MSE(alpha_CWT*CWT,WVS);
MSE_vec_CWT = squeeze(MSE_vec_CWT);
mMSE_CWT = mean(MSE_vec_CWT); 
std_MSE_CWT = std(MSE_vec_CWT);


%%  Plot an example using optimal spectral parameters for each method

r0 = randi([1,num_real]); % random realization for plotting an example
% r0=11;

exampleS_HANN = S_HANN(:,:,r0,I_HANN_opt); 
exampleS_WOSA = S_WOSA(:,:,r0,I_WOSA_opt); 
exampleS_LSP_HATS = S_LSP_HATS(:,:,r0); 
exampleS_LSP_true_par = S_LSP_true_par(:,:,r0); 

figure

subplot(321)
c=[min(min(WVSshift)) max(max(WVSshift))];
pcolor(TI,FI,WVSshift')  
shading interp
caxis(c)
ylim ([0 100])
ylabel('Frequency (Hz)')
xlabel('Time (s)')
title(sprintf('WVS with f_0= %d Hz', f0));
colorbar

subplot(322)
c=[min(min(exampleS_HANN)) max(max(exampleS_HANN))];
pcolor(TI,FI,exampleS_HANN')  
ylim ([0 100])
shading interp
caxis(c)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar
title(sprintf('HANN (win length = %d)',win_HANN_opt))

subplot(323)
c=[min(min(exampleS_WOSA)) max(max(exampleS_WOSA))];
pcolor(TI,FI,exampleS_WOSA')  
ylim ([0 100])
shading interp
caxis(c)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar
title(sprintf('WOSA (K = %d)',win_WOSA_opt))

subplot(324)
c=[min(min(exampleS_LSP_HATS)) max(max(exampleS_LSP_HATS))];
pcolor(TI,FI,exampleS_LSP_HATS')  
ylim ([0 100])
shading interp
caxis(c)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar
title('LSP HATS')

subplot(325)
c=[min(min(exampleS_LSP_true_par)) max(max(exampleS_LSP_true_par))];
pcolor(TI,FI,exampleS_LSP_true_par')  
ylim ([0 100])
shading interp
caxis(c)
ylabel('Frequency (Hz)')
xlabel('Time (s)')
colorbar
title('LSP true par')

sgtitle('Example')

%% Boxplots

if (num_real>1)   

    s1 = mMSE_HANN_matrix(I_HANN_opt,:); 
    s2 = mMSE_WOSA_matrix(I_WOSA_opt,:);
    s3 = MSE_vec_WV;
    s4 = reshape(MSE_LSP_true_par,1,[]);
    s5 = reshape(MSE_LSP_HATS,1,[]); % reshape to have the values in a vector
    s6 = MSE_vec_CWT; 

    
    figure
    h = boxplot([s1' s2' s3 s4' s5' s6],'labels',{'HANN','WOSA','WV','LSP','LSP-HATS','CWT'});
    title('Boxplot of the MSE')
    set(h,'linewidth',1.3) 
    ylabel('MSE')

end

disp(sprintf('MSE HANN %f (window length = %d)', mMSE_HANN_opt, win_HANN_opt ))
disp(sprintf('MSE WOSA %f (%d windows)', mMSE_WOSA_opt, win_WOSA_opt))
disp(sprintf('MSE WV %f', mMSE_WV))
disp(sprintf('MSE LSP true parameters %f', mMSE_LSP_true_par))
disp(sprintf('MSE LSP-HATS %f', mMSE_LSP_HATS))
disp(sprintf('MSE CWT %f', mMSE_CWT))
