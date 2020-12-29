function [X,X_freq,C,C_freq,R,R_freq,Q]  =  lsp_f0_sim(num_sim, f0 ,a_q, b_q, c_q, c_r,sigma_square,N,time);

% lsp_f0_sim computes num_sim data-vector of a locally stationary process 
% with the model parameters specified and centre frequency f0
%   
% OUTPUT:
%   X: output data vector [N x num_sim]
%   X(:,i) represent the i-th realization
%   C: Covariance matrix
%   R: Toeplitz matrix (stationary part of covariance)
%   Q: Hankel matrix (time varying part of covariance)

% INPUT:
%   num_sim: number of simulated trajectories
%   f0: centre frequency
%   a_q, b_q, c_q, c_r: Model parameters
%   sigma_square: constant noise level
%   N: one trajectory length
%   time: time vector of length N
%

t = time;
s = t';

tau_R = t*ones(1,length(s))-ones(length(t),1)*s;
tau_Q = (t*ones(1,length(s))+ones(length(t),1)*s)./2;

R = exp(-(c_r/8).*(tau_R).^2);
R_freq = exp(-(c_r/8).*(tau_R).^2).*cos(2*pi*f0*tau_R);
Q = sigma_square+a_q*exp(-(c_q/2)*((tau_Q-b_q).^2));

% Covariance matrix of the base-band lsp-process
C = (R.*Q);
C_freq = (R_freq.*Q);

% Realizations from the filtered white noise realization b
Noise = randn(N,num_sim); % Noise realization
c1 = sqrtm(C);
c2 = sqrtm(C_freq);
X  = real(c1)*Noise;
X_freq = real(c2)*Noise; %lsp-realization



    