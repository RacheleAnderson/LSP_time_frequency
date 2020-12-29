function [noise_HATS, a_q_HATS, b_q_HATS, c_q_HATS, c_r_HATS]= HATS(X, time, LB, UB, THETA_0);

% Inference using HATS on a LSP defined through a covariance function
% C(s,t)= q((s+t)/2)*r(t-s) with q,r functions:
% q(x) = L + a_q* exp(-(c_q/2) * ((x-b_q)^2))
% r(x) =  exp(-(c_r/8).*x.^2)

% INPUT:
%   X: matrix of size N x numReal containing process 'numReal' realizations each of length N
%   time: time vector
%   LB, UB, THETA_0 : lower bounds, upper bounds, starting values for optimization

% OUTPUT:
%   L_est, a_q_est, b_q_est, c_q_est, c_r_est : estimated LSP parameters
    
N = size(X,1);
M = size(X,2); %number of realizations

T0 = time(1);
Tf = time(end);
t = time;
s = t';
tau_Q = (t*ones(1,length(s))+ones(length(t),1)*s)./2; % matrix
tau_R = t*ones(1,length(s))-ones(length(t),1)*s;
tauQ = diag(tau_Q); % vector
tauR = diag(flipud(tau_R));

% optimization settings from input:
noise_LB = LB(1); a_q_LB = LB(2); b_q_LB = LB(3); c_q_LB = LB(4); c_r_LB = LB(5);
noise_UB = UB(1); a_q_UB = UB(2); b_q_UB = UB(3); c_q_UB = UB(4); c_r_UB = UB(5);
noise_0 = THETA_0(1); a_q_0 = THETA_0(2); b_q_0 = THETA_0(3); c_q_0 = THETA_0(4); c_r_0 = THETA_0(5);

% SCM
if (M > 1)
    C_SC = cov(X');
else
    C_SC = X*X';
end
% PART I: estimate Q

all_AIP = X.^2; % all_AIP (:,i) contains the average istanteneous power of realization i

meanAIP = mean(all_AIP,2); % Average istantaneous power for the process 

% Estimate the lambda parameters by non-linear LS fitting of the function q
qfo = fitoptions('Method','NonlinearLeastSquares','Lower',[noise_LB a_q_LB b_q_LB c_q_LB],'Upper',[noise_UB a_q_UB b_q_UB c_q_UB],'StartPoint',[noise_0 a_q_0 b_q_0 c_q_0]);
qft = fittype(' L + a* exp(-(c/2) * ((x-b).^2))','options',qfo);
q_fit = fit(tauQ,meanAIP,qft); 

% Round off the estimates to the second decimal
noise_HATS = round (q_fit.L * 100) / 100;
a_q_HATS = round(q_fit.a * 100) / 100;
b_q_HATS = round(q_fit.b * 100) / 100;
c_q_HATS = round(q_fit.c * 100) / 100;

Q_HATS = noise_HATS + a_q_HATS * exp(-(c_q_HATS/2)*((tau_Q-b_q_HATS).^2));  % estimated matrix Q

% PART II: Estimation of stationary covariance function r
R_0 = zeros(N,N);
R_SC = C_SC ./ (Q_HATS);
for k=1: N % we have N diagonals and the index k represents how long is the diagonal

    pos = N-k; % position of the diagonal (i.e. shift from the main diagonal)
    mean_diag = mean (diag(R_SC, pos));     
    aux = diag (mean_diag * ones(1,k), pos); 
    R_0 = R_0 + aux;

    clearvars pos mean_diag aux
end
% Make it symmetric
aux = R_0';
aux(logical(eye(size(aux)))) = 0;
R_0 = R_0 + aux;

r_antidiag_HATS = diag(flipud(R_0));
R_fit_HATS = r_antidiag_HATS ;

% Estimate the rho parameters by LS fitting of the function r
rfo = fitoptions('Method','NonlinearLeastSquares','Lower',[c_q_HATS],'Upper',[c_r_UB],'StartPoint',[c_r_0]);
rft = fittype('exp(-(c/8).*x.^2)','options',rfo);
r_HATS = fit(tauR,R_fit_HATS,rft);
c_r_HATS = round(r_HATS.c * 100) / 100; % round off the estimates to the second decimal

% Generate the matrix R and final covariance estimate
R_HATS = exp(-(c_r_HATS/8).*(tau_R).^2);
C_tot_HATS = (R_HATS.*Q_HATS);

% Check if the matrix is positive semi-definite
A = C_tot_HATS;
eig_A = eig(A);
eig_A_round = round(eig_A * 1000)/1000;

if (min(eig_A_round)<0)
   display('Error: covariance estimate is not positive semi-definite!')
end


end
