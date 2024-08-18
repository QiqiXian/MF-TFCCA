function result = VAR_est1(Data, p, h, lambda)
% PURPOSE: Fit (p,h)-autoregression (i.e. VAR(p) model iterated h-times) 
%          with Newey and West's (1987) HAC estimator. Newey and West's (1994) 
%          automatic bandwidth selection is available.
%---------------------------------------------------------------------
% USAGE: result = VAR_est1(Data, p, h, lambda)
% where: Data = data (T-by-K-by-N matrix)
%        T = # of observations
%        K = dimension of VAR
%        N = # of trials
%        p = # of lags            (default = 1)
%        h = prediction horizon   (default = 1) 
%        lambda = tuning parameter for Newey-West HAC
%                 (default = 'NW' i.e. Newey and West's automatic lag selection)
%----------------------------------------------------------------------
% RETURNS: result.Data    = Data
%          result.OLSE    = pK-by-K matrix containing OLSE
%          result.Omegah  = K-by-K covariance matrix for error
%                           in (p,h)-autoregression
%          result.se      = element-wise standard errors for OLSE  (pK-by-K)
%          result.Sigma   = Newey-West HAC estimator for OLSE cov. mat. (pK^2-by-pK^2)
%          result.horizon = h
%          result.nlag    = p
%          result.nobs    = T
%          result.ndim    = K
%          result.T_star  = T - p - h + 1 (effective sample size)
%          result.residmat = residual in matrix form  (T_star-by-K)
%          result.pdflag   = 1 if Delta_hat_p_h is positive definite
%                            0 otherwise
%          result.n        = bandwidth
%          result.OLSE1    = pK-by-K matrix containing OLSE for h = 1
%                            (necessary for bootstrapped causality test and impulse response)
%          result.Omega1   = K-by-K covariance matrix for error
%                            in (p,1)-autoregression
%                            (necessary for bootstrapped causality test and impulse response)
%          result.lambda   = lambda
%          result.analcov_flag = 1 if analytical covariance is substituted
%                                     to Sigma
%                                0 otherwise
%                                (default: 0. This equals 1 only if analytical
%                                 covariance is used in benchmark simulation.)
%          result.stab         = 1 if stability condition for VAR(p) holds.
%                                0 otherwise.
% ----------------------------------------------------------------------
% Note: Be careful about the order of elements in OLSE. 
%       k-th block of OLSE contains the TRANSPOSE of A_k(h) for
%       k = 1,...,p
%       (p,h)-autoregression:
%       X(t) = sum_{k=1}^{p} A_k(h) X(t-h-k+1) + epsilon(t).
% 
%       This code assumes no constant terms. Demean each series in advance.
% ----------------------------------------------------------------------
% References: E. Ghysels, J. B. Hill, & K. Motegi (2013)
%             Testing for Granger Causality with Mixed Frequency Data. 
%             Working Paper at UNC. 
% 
%             J. M. Dufour, D. Pelletier, & E. Renault (2006)
%             Short Run and Long Run Causality in Time Series: Inference.
%             Journal of Econometrics, 132, 337-362.
%
%             W. K. Newey & K. West (1987)
%             A Simple, Positive Semi-definite, Heteroskedasticity and 
%             Autocorrelation Consistent Covariance Matrix.
%             Econometrica, 55, 703-708.
%
%             W. K. Newey & K. West (1994)
%             Automatic Lag Selection in Covariance Matrix Estimation.
%             Review of Economic Studies, 61, 631-653.
%       
%             J. D. Hamilton (1994) Time Series Analysis.
%             New Jersey: Princeton University Press.
% -------------------------------------------------------------------
% Written by Kaiji Motegi, Waseda University.
% Sep. 9, 2014.
% ----------------------------------------------------------------------

% K = dimension of VAR
% T = # of observations
% N = # of trials
Data = permute(Data, [2,1,3]);
[K, T, ntrials] = size(Data);

% Default settings
if nargin == 1
    p = 1; h = 1; lambda = 'NW';
elseif nargin == 2
    h = 1; lambda = 'NW';
elseif nargin == 3
    lambda = 'NW';
end;

% OLSE for h = 1 will be needed for Bootstrap whatever actual h is.
if h > 1
    % construct dependent variable
%     W_1_1 = Data((p+1):T, :);
    W_1_1 = reshape(Data(:, (p+1):T, :), K, (T-p)*ntrials)';

    % construct independent variables (i.e. lagged dependent variables)
    % [x1(-1), ..., xK(-1), x1(-2), ..., xK(-2), ..., x1(-p), ..., xK(-p)]
    % no constants. center each variable before using this code.
%     W_bar_p_1 = zeros(T-p, p * K);
    W_bar_p_1 = zeros((T-p)*ntrials, p * K);
    for i = 1:p     % lag i
         first = (i-1)*K + 1;
         last = i*K;
%          W_bar_p_1(:, first:last) = Data((p+1-i):(T-i), :);
         W_bar_p_1(:, first:last) = reshape(Data(:,(p+1-i):(T-i),:), K, (T-p)*ntrials)';
    end;    

    % get OLS estimator (pK x K)
    B_hat_1 = W_bar_p_1 \ W_1_1;
    
    % residual matrix ((T-p) x K)
    U_1_hat_1 = W_1_1 - W_bar_p_1 * B_hat_1;

    % compute covariance matrix estimator of error 
    Omega_1_hat = cov(U_1_hat_1);
end;    

% Effective sample size
T_star = T - p - h + 1;

% Construct dependent variables in matrix form
% W_h_h = Data((p+h):T, :);
W_h_h = reshape(Data(:, (p+h):T, :), K, T_star*ntrials)';

% Construct independent variables (i.e. lagged dependent variables)
% [x1(-h), ..., xK(-h), x1(-h-1), ..., xK(-h-1), ..., x1(-h-p+1), ..., xK(-h-p+1)]
% no constants. center each variable before using this code.
% W_bar_p_h = zeros(T_star, p*K);
W_bar_p_h = zeros(T_star*ntrials, p*K);
for i = 1:p
     first = (i-1)*K + 1;
     last = i*K;
%      W_bar_p_h(:, first:last) = Data((p-i+1):(T_star+p-i), :);
     W_bar_p_h(:, first:last) = reshape(Data(:,(p-i+1):(T_star+p-i),:), K, T_star*ntrials)';
end;    

% Get OLS estimator
B_hat_h = W_bar_p_h \ W_h_h;

% residual matrix
U_h_hat_h = W_h_h - W_bar_p_h * B_hat_h;

% covariance matrix of redidual  
Omega_h_hat = cov(U_h_hat_h);

% Construct some matrices for covariance estimator

% covariance matrix of independent variables
Gamma_hat_p_0 = cov(W_bar_p_h);

% independent variables times residuals
Y_hat_p = zeros(p*K^2, T_star);
for t = 1:T_star
     Y_hat_p(:,t) = reshape(W_bar_p_h(t,:)' * U_h_hat_h(t,:), p*K^2, 1);
end;

% OPTIONAL: choose automatic tuning parameter based on Newey and West (1994)
if strcmp(lambda,'NW')   
    x = Y_hat_p(:, 1:(end-1));
    y = Y_hat_p(:, 2:end);
    A = x \ y;
    H = y - x * A;
    N = floor(4 * (T_star / 100)^(2/9));
    temp = (ones(1, p*K^2) * H)';
    sigma = zeros(N+1, 1);
    for j = 1:(N+1)
         sigma(j) = 1/(T_star - 1) * sum(temp(j:(end-1)).*temp(1:(end-j)));
    end;
    s1 = 2 * (1:N) * sigma(2:end);
    s0 = sigma(1) + 2 * sum(sigma(2:end));
    lambda = 1.1447*((s1/s0)^2)^(1/3);
end;

% bandwidth
n = min(max(h, floor(lambda * T_star^(1/3))), T_star);

% Compute long-run variance of Y based on Bartlett kernel
Delta_hat_p_0_h = zeros(p*K^2, p*K^2);
for t = 1:T_star
    Delta_hat_p_0_h = Delta_hat_p_0_h + Y_hat_p(:,t) * Y_hat_p(:,t)';
end;
Delta_hat_p_0_h = Delta_hat_p_0_h / T_star;

Delta_hat_p_h = Delta_hat_p_0_h;
if n > 1
   % s-th order autocovariance of Y for s = 1, ..., n-1
   % The case that s = n does not have to be considered since it has zero weights.
   Delta_hat_p_s_h = zeros(p*K^2, p*K^2, n);
   for s = 1:(n-1)
       for t = (s+1):T_star
            Delta_hat_p_s_h(:,:,s) = Delta_hat_p_s_h(:,:,s) + Y_hat_p(:,t) * Y_hat_p(:,t-s)';
       end;
       Delta_hat_p_s_h(:,:,s) = Delta_hat_p_s_h(:,:,s) / T_star;
   end;

   % vector of weights (Bartlett)
   w = 1 - (1:(n-1))' / n;

   % Add weighted matrices
   for s = 1:(n-1)
        Delta_hat_p_h = Delta_hat_p_h + w(s)*(Delta_hat_p_s_h(:,:,s) + Delta_hat_p_s_h(:,:,s)');
   end;
end;    

% pdflag equals 1 if the long-run variance of Y is positive definite and 0 o.w.
pdflag = all(eig(Delta_hat_p_h) > 0);

% Newey and West's (1987) HAC estimator for OLSE covariance matrix
temp = kron(eye(K), inv(Gamma_hat_p_0));
Sigma_hat_p_B_hat_h = temp * Delta_hat_p_h * temp;

% element-wise standard errors for OLSE
semat = reshape(sqrt(diag(Sigma_hat_p_B_hat_h) / T), p*K, K);

% save results
result.Data = Data;
result.OLSE = B_hat_h;
result.Omegah = Omega_h_hat;
result.se = semat;
result.Sigma = Sigma_hat_p_B_hat_h;
result.horizon = h;
result.nlag = p;
result.nobs = T;
result.ndim = K;
% result.T_star = T_star;
result.T_star = T_star * ntrials;
result.residmat = U_h_hat_h;
result.pdflag = pdflag;
result.n = n;
result.lambda = lambda;
result.analcov_flag = 0;
if h == 1
    result.OLSE1 = B_hat_h;
    result.Omega1 = Omega_h_hat;
else
    result.OLSE1 = B_hat_1;
    result.Omega1 = Omega_1_hat;
end;    

% check stability condition
OLSE1 = result.OLSE1;
% construct F: see Eq. (10.1.10) in p.259, Hamilton (1994). 
F = zeros(p*K, p*K);
for i = 1:p
    first = (i-1) * K + 1;
    last = i * K;
    F(1:K, first:last) = OLSE1(first:last, :)';
end;    
if p > 1
    F((K+1):(p*K), 1:((p-1)*K)) = eye((p-1)*K, (p-1)*K);
end;    
% If the eigenvalues of F all lie inside the unit circle, VAR is stationary
result.stab = all(abs(eig(F)) < 1);