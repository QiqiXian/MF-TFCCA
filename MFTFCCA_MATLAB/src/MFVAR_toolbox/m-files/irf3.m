function [IRF, lb, ub] = irf3(result_VAR, hmax, figureflag, alpha, bsnum, labels)
% PURPOSE: Compute impulse response function at horizon 0, 1, ..., hmax
%          along with bootstrapped confidence intervals. 
%          Cholesky decomposition is used.
%---------------------------------------------------------------------
% USAGE: [IRF, lb, ub] = irf3(result_VAR, hmax, figureflag, alpha, bsnum, labels)
% where: result_VAR = structure from VAR_est.m
%              hmax = maximum horizon  (scalar)
%        figureflag = 1 if figure is desired.
%                     0 otherwise.
%                     (default = 1)
%             alpha = nominal size  (scalar)
%                     100(1-alpha)% condifence interval is computed via bootstrap.
%                     (default = 0.05)
%             bsnum = # of bootstrap replications  (scalar)
%            labels = vector of name labels for a figure (K x 1)
%                     (default = [var1; var2;, ...])
%                 K = dimension of VAR
%----------------------------------------------------------------------
% RETURNS: IRF = impulse response functions (K x K x (hmax+1))
%                IRF(i,j,h) represents (h-1)-step-ahead response of variable i to 
%                1-sigma shock of variable j, where h = 1, ..., hmax + 1.
%          lb  = lower bounds of confidence intervals (K x K x (hmax+1))
%                lb(i,j,h) represents the lower bound of IRF(i,j,h).
%          ub  = upper bounds of confidence intervals (K x K x (hmax+1))
%                ub(i,j,h) represents the upper bound of IRF(i,j,h).
% ----------------------------------------------------------------------
% References: J. Hamilton (1994). Time Series Analysis. Princeton
%                University Press. Section 11.4.
%             H. Lutkepohl (2005). New Introduction to Multiple Time Series
%                Analysis. Springer-Verlag, Berlin. Section 3.7.
%             H. Lutkepohl, A. Staszewska-Bystrova, and P. Winker (2013).
%                Comparison of Methods for Constructing Joint Confidence
%                Bands for Impulse Response Functions. MAGKS Joint
%                Discussion Paper Series in Economics.
% -------------------------------------------------------------------
% Written by Kaiji Motegi, Waseda University.
% Aug. 8, 2014.
% ----------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sub-function that computes impulse response (not confidence intervals) 
function [IRF, Amat, Omega] = impulse(result_VAR, hmax)
          % covariance matrix of error term of (p,1)-autoregression 
          Omega = result_VAR.Omega1;
          % Cholesky decomposition. P is a lower triangular such that Omega = PP'.
          P = chol(Omega, 'lower');

          % # of lags
          nlag = result_VAR.nlag;
          % dimension of VAR
          ndim = result_VAR.ndim;

          % construct a ndim x ndim x nlag matrix of OLS estimates from result_VAR.OLSE1.
          % Note how result_VAR.OLSE1 is constructed.
          % It is a (nlag * ndim) x ndim matrix containing least squares estimates for
          % (nlag,1)-autoregression.
          % The first ndim x ndim block is the TRANSPOSE of A_1 and so on.
          % Amat removes the transpose.
          Amat = zeros(ndim, ndim, nlag);
          for k = 1:nlag
               first = (k-1)*ndim + 1;
               last = k*ndim;
               Amat(:,:,k) = result_VAR.OLSE1(first:last, :)';
          end;    

          % Construct VMA(infty) coefficients up to lag hmax
          % Recursive formula: Psi(k) = sum_{s=1}^{nlag} A(s) Psi(k-s)         
          Psimat = zeros(ndim, ndim, hmax);
          indmat = repmat((1:hmax)', 1, nlag) - repmat(1:nlag, hmax, 1);  % index matrix
          % indmat = [   1 - 1, ... ,     1 - nlag;
          %                :  ,  :  ,       :     ;
          %           hmax - 1, ... ,  hmax - nlag];
          
          for k = 1:hmax   % compute Psi(k)
               for l = 1:nlag      % add nlag terms
                    ind = indmat(k,l);
                    if ind < 0     % Psi with negative index is a null matrix
                        increment = zeros(ndim, ndim);
                    elseif ind == 0   % Psi with zero index is eye(K)
                        increment = Amat(:,:,l);
                    else           % Psi with positive index is just as is
                        increment = Amat(:,:,l) * Psimat(:,:,ind);
                    end;    
                    Psimat(:,:,k) = Psimat(:,:,k) + increment;
               end;     
          end;    

          % compute impulse response functions
          % IRF(i,j,h) represents h-step-ahead response of variable i to 
          % 1-sigma shock of variable j. 
          IRF = zeros(ndim, ndim, hmax+1);
          IRF(:,:,1) = P;  % 0-step ahead IRF
          for h = 1:hmax
               IRF(:, :, h+1) = Psimat(:,:,h) * P;
          end;    
end  %   close sub-function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Main function that computes confidence band and draws a figure
% default setting 
if nargin == 2
    figureflag = 1;   % draw a figure
    alpha = 0.05;     % construct 95% confidence interval
    bsnum = 1000;     % 1000 bootstrap replications for each IRF
elseif nargin == 3
    alpha = 0.05;
    bsnum = 1000;
elseif nargin == 4
    bsnum = 1000;
end;    

% get impulse response and least squares estimates from sub-function
[IRF, Amat, Omega] = impulse(result_VAR, hmax);
p = result_VAR.nlag;          % VAR lag order
K = result_VAR.ndim;          % VAR dimension
T = result_VAR.nobs;          % sample size
lambda = result_VAR.lambda;   % tuning parameter for Newey-West HAC estimator

% begin bootstrap to get artificial IRFs
IRF_bs = zeros(K, K, hmax+1, bsnum);
for b = 1:bsnum
     % simulate VAR processes using OLS estimates
     E = mvnrnd(zeros(1,K), Omega, T);
     X_bs = sim_VAR(E, Amat);    
     
     % fit VAR(p) and get IRF
     result_VAR_bs = VAR_est1(X_bs, p, 1, lambda);
     [IRF_bs(:,:,:,b), ~, ~] = impulse(result_VAR_bs, hmax);
end;    

% sort IRF, fixing variables and horizons
IRF_bs_sort = sort(IRF_bs, 4);

lb = IRF_bs_sort(:,:,:, floor(0.5*alpha*bsnum));         % lower bound
ub = IRF_bs_sort(:,:,:, floor((1 - 0.5*alpha)*bsnum));   % upper bound

if figureflag == 1  % draw figure
    figure
    for i = 1:K    % response of variable i
         for j = 1:K       % 1-sigma shock in variable j
              % call IRF, lower bound, and upper bound
              IRF_ij = reshape(IRF(i,j,:), hmax+1, 1);
              lb_ij = reshape(lb(i,j,:), hmax+1, 1);
              ub_ij = reshape(ub(i,j,:), hmax+1, 1);
              
              % draw a panel.
              % blue, solid line for IRF.
              % red, dashed lines for confidence intervals.
              index = (i-1)*K + j;              
              subplot(K, K, index)
              linescale = 0:hmax; % 0-step ahead through hmax-step ahead IRF
              plot(linescale, IRF_ij, 'b-', linescale, lb_ij, 'r--', linescale, ub_ij, 'r--');
              xlim([0, hmax]);
              if i == K 
                  if nargin == 6   % specific names
                      % strcat deletes empty spaces
                      xlabel(['1\sigma shock in ', strcat(labels(j,:))]);
                  else             % default names
                      xlabel(['1\sigma shock in var.', num2str(j)]);
                  end;    
              end;    
              if j == 1
                  if nargin == 6   % specific names
                      ylabel(['IR of ', strcat(labels(i,:))]);
                  else             % default names
                      ylabel(['IR of var.', num2str(i)]);
                  end;    
              end;    
              hold on
         end;
    end;     
    hold off
end;    

end   % close the main function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%