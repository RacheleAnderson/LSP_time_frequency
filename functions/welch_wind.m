function [wind, wei]=welch_wind(N,K);

% Calculates the 'now' (default 1) windows 
% 
% K : A number between 1 and 8 (default 1)
% N   : Window length
% wind: The windows in a matrix as columns 
% wei: Vector of weights

  NN = ceil(2/(K+1)*N); % actual window length
  step=round(1/(K+1)*N); 
  win = hanning(NN);
  win = win./sqrt(win'*win);
  wei = ones(K,1)/K; %1/sqrt(K);
  wind=[win; zeros(N-NN,1)];
  if K>=2
    for i=2:K-1
      wind=[wind [zeros(step*(i-1),1) ;win; zeros(N-NN-step*(i-1),1)]];
    end
    wind=[wind [zeros(N-NN,1); win]];
  end  
  
end
