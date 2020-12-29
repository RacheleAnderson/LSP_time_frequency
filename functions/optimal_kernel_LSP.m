function [uopt,sopt,fiopt,FFiopt]  =  optimal_kernel_LSP(NN, a_q, c_q, c_r, L, fs);

% optimal_kernel_LSP computes the mse-optimal multitapers and weights corresponding
% to a LSP with covariance defined as in 
% Anderson, R., Sandsten, M. "Time-frequency feature extraction for classification of episodic memory"
% EURASIP J. Adv. Signal Process. 2020, 19 (2020). https://doi.org/10.1186/s13634-020-00681-8
%
% LSP model:
%   q(x)  =  L + a_q* exp(-(c_q/2) * ((x-b_q)^2))
%   r(x)  =   exp(-(c_r/8).*x.^2)
% Note: the parameter b_q is not relevant in time-frequency analysis since it does not affect the kernel

% OUTPUT:
%   uopt: LSP-optimal multitapers
%   sopt: LSP-optimal weights
%   fiopt: LSP-optimal ambiguity kernel
%   FFiopt: LSP-optimal time-frequency kernel
%
% INPUT:
%   NN: multitaper window lengths
%   a_q, c_q, c_r, L : parameters of the process (L  =  noise level)


LL = 2*NN;
N = NN/2;

tl = [-NN:NN-1]'/fs;
f = fs*[-NN:NN-1]'/(LL);

% The optimal LSP-kernel is:
% (1/(2*pi))*[ abs(Q(theta))^2 * abs(r(tau))^2 ] / [abs(Q(theta))^2 * abs(r(tau))^2 + (FT(abs(r))^2)(theta) * (IF(abs(Q))^2(tau))]
%  =  1 / (1 + STAR) where
% STAR  =  [(FT(abs(r))^2)(theta) * (IF(abs(Q))^2(tau))] / [ [abs(Q(theta))^2 * abs(r(tau))^2 ]

theta = 2*pi*f*ones(1,length(tl));
tau = ones(length(f),1)*(tl');

r_square = exp(-(c_r/4)*tau.^2); % in tau
ft_r_square = 2*sqrt(pi/c_r)*exp(-(1/c_r)*theta.^2);% in theta

discreteDirac = dirac(theta);
discreteDirac(isinf(discreteDirac)) = 1;
Q_square = L^2*discreteDirac+2*pi*a_q^2/c_q*exp(-(1/c_q)*theta.^2) + 2*a_q*sqrt(2*pi/c_q)*L*discreteDirac.*exp(-(1/(2*c_q))*theta.^2);

ift_Q_square = (1/(2*pi))*(L^2 + 2*pi*sqrt(pi/c_q)*a_q^2*exp(-(c_q/4)*tau.^2)+ 2*pi*L*a_q*sqrt(2*pi/c_q)); % in tau

Qr2 = Q_square.*r_square; 
Fr2iQ2 = ft_r_square.*ift_Q_square;

fiopt = 1./ (1+ (Fr2iQ2./Qr2));

indexofNaN = find(isnan(fiopt));  % if we are exceeding matlab precision 
fiopt(indexofNaN) = 0;           % but we know that it should be 0

Fiopt = zeros(LL,LL);
Giopt = zeros(LL,NN);

% Optimal time-lag-kernel
for i = 1:LL
    Giopt(:,i) = fftshift((ifft([fiopt(NN+1:LL,i);fiopt(1:NN,i)])));
end
Giopt = Giopt(NN-N+1:NN+N,:);

% Optimal doppler-kernel
for i = 1:LL
    Fiopt(i,:) = fftshift((fft([fiopt(i,NN+1:2*NN) zeros(1,LL-2*NN) fiopt(i,1:NN)]))); 
end

% Optimal TF-kernel
for i = 1:LL
    FFiopt(:,i) = real(fftshift((ifft([Fiopt(NN+1:LL,i);Fiopt(1:NN,i)]))));
end


t = [-N:N-1]';

% Optimal multitapers uopt and weights sopt computed from the rotated time-lagkernel

% NB. The multiple windows are obtained as the eigenvectors of the rotated
% time-lag estimation kernel. The spectrograms from the different windows
% are weighted with the eigenvalues and the resulting multiple window
% spectrogram is an estimate of the optimal smoothed WVS

for i = 1:NN
    for ii = 1:NN
        if abs((i-1)-(ii-1))/2-fix(abs((i-1)-(ii-1))/2)<0.1
            
            RRR(ii,i) = Giopt(((i-1)+(ii-1))/2+1,NN+(i-1)-(ii-1)+1);
            
        end    
     end
end

for i = 1:NN
    for ii = 1:NN
        if abs((i-1)-(ii-1))/2-fix(abs((i-1)-(ii-1))/2)>0.1
            if i>1 & ii>1 & i<NN & ii<NN
               RRR(i,ii) = (RRR(i-1,ii)+RRR(i+1,ii)+RRR(i,ii-1)+RRR(i,ii+1))/4;
            elseif i==1 & ii>1 & ii<NN
               RRR(1,ii) = (RRR(i+1,ii)+RRR(1,ii-1)+RRR(1,ii+1))/3;
            elseif ii==1 & i>1 & i<NN
               RRR(i,1) = (RRR(i-1,1)+RRR(i+1,1)+RRR(i,ii+1))/3;
            elseif i==NN & ii>1 & ii<NN
               RRR(NN,ii) = (RRR(i-1,ii)+RRR(NN,ii-1)+RRR(NN,ii+1))/3;
            elseif ii==NN & i>1 & i<NN
               RRR(i,NN) = (RRR(i-1,NN)+RRR(i+1,NN)+RRR(i,ii-1))/3;
            end
        end    
     end
end

[u1,s] = eig(RRR);
s = diag(real(s));
s_max = max(s);
k = find(s>0.1*s_max); 
uopt = u1(:,k);
sopt = s(k);


