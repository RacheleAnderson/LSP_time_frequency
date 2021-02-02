
function MSE = compute_MSE(S_est,true_WVS)

% Input:
% S_est: spectrum estimate
% true_WVS: true Wigner-Ville Spectrum

% Output:
% MSE: mean square error 

SE = abs(S_est-true_WVS).^2;
MSE_f = mean(SE,2); % mean on freq
MSE = mean(MSE_f,1); % mean on time

end