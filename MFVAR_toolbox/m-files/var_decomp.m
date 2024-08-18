function outcome = var_decomp(result_VAR, hmax)
% PURPOSE: Forecast error variance decomposition at horizon 0, ..., hmax-1.
%          Cholesky decomposition is used.
%---------------------------------------------------------------------
% USAGE: outcome = var_decomp(result_VAR, hmax)
% where: result_VAR = structure from VAR_est.m
%        hmax = maximum horizon
%----------------------------------------------------------------------
% RETURNS: outcome = variance decomposition (K x K x hmax)
%                     where K = VAR dimension
%          outcome(j,k,h) represents how much portion variable k
%          contributes to the total variation of variable j at horizon h-1.
%          sum_{k = 1}^{K} outcome(j,k,h) = 1 for any j and h.
% ----------------------------------------------------------------------
% References: H. Lutkepohl (2005). New Introduction to Multiple Time Series
%                Analysis. Springer-Verlag, Berlin. Section 2.3.3.
% -------------------------------------------------------------------
% Written by Kaiji Motegi, Waseda University.
% Aug. 9, 2014.
% ----------------------------------------------------------------------

% covariance matrix of error term of (p,1)-autoregression 
Omega = result_VAR.Omega1;
% Cholesky decomposition. P is a lower triangular such that Omega = PP'.
P = chol(Omega, 'lower');

% # of lags
p = result_VAR.nlag;
% dimension of VAR
K = result_VAR.ndim;

% construct a K x K x p matrix of OLS estimates from result_VAR.OLSE1.
% Note how result_VAR.OLSE1 is constructed.
% It is a pK x K matrix containing least squares estimates for
% (p,1)-autoregression.
% The first K x K block is the TRANSPOSE of A_1 and so on.
% Amat removes the transpose.
Amat = zeros(K, K, p);
for k = 1:p
     first = (k-1)*K + 1;
     last = k*K;
     Amat(:,:,k) = result_VAR.OLSE1(first:last, :)';
end;    

% Construct VMA(infty) coefficients up to lag hmax
% Recursive formula: Psi(k) = sum_{s=1}^{p} A(s) Psi(k-s)
Psimat = zeros(K, K, hmax);
indmat = repmat((1:hmax)', 1, p) - repmat(1:p, hmax, 1);  % index matrix
for k = 1:hmax   % compute Psi(k)
     for l = 1:p      % add p terms
          ind = indmat(k,l);
          if ind < 0     % Psi with negative index is a null matrix
              increment = zeros(K,K);
          elseif ind == 0   % Psi with zero index is eye(K)
              increment = Amat(:,:,l);
          else           % Psi with positive index is just as is
              increment = Amat(:,:,l) * Psimat(:,:,ind);
          end;    
          Psimat(:,:,k) = Psimat(:,:,k) + increment;
     end;     
end;    

% multiply VMA coefficients and Cholesky factor
Thetamat = zeros(K, K, hmax);
for k = 1:hmax
     Thetamat(:,:,k) = Psimat(:,:,k) * P;
end;    

% variance decomposition
I = eye(K);
contribution = zeros(K, K, hmax);
outcome = zeros(K, K, hmax);
for h = 1:hmax   % fix prediction horizon (strictly, horizon under consideration is h-1)
     for k = 1:K   % consider the contribution of variable k to variable j
          for j = 1:K   
               % See Eq. (2.3.36) in Lutkepohl (2005).
               contribution(j,k,h) = (I(:,j)' * P * I(:,k))^2;
               if h > 1
                   for i = 1:(h-1)
                        contribution(j,k,h) = contribution(j,k,h) + ...
                                               (I(:,j)'*Thetamat(:,:,i)*I(:,k))^2;
                   end;     
               end;    
          end;   
     end;   
     % total variation of variable j at horizon h-1
     MSE = sum(contribution(:,:,h), 2);
     % calculate portions
     outcome(:,:,h) = contribution(:,:,h) ./ repmat(MSE, 1, K);
end;    
