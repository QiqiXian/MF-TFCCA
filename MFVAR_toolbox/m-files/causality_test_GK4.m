function [W, P_value] = causality_test_GK4(result_VAR_est, from, to, iota, iotap, bs)
% PURPOSE: Test for Granger noncausality from "from" variables to 
%          "to" variables in VAR models
%--------------------------------------------------------------------
% USAGE: [W, P_value] = causality_test_GK3(result_VAR_est, from, to, iota, iotap, bs)
% where: result_VAR_est = result of VAR_est.m  (structure)
%                  from = index for "from" variables 
%                         (scalar or 2-dimensional vector)
%                    to = index for "to" variables
%                         (scalar or 2-dimensional vector)
%                  iota = increment for "to" variables
%                         (scalar; default = 1)
%                 iotap = increment for "from" variables
%                         (scalar; default = 1)
%                    bs = # of bootstrap replications
%                         'off' if you do not use bootstrap
%                         (default = 'off')
%--------------------------------------------------------------------
% RETURNS: W       = Wald statistic
%          P_value = Associated P value 
%--------------------------------------------------------------------
% References: E.Ghysels, J.B.Hill, & K.Motegi (2013)
%             Testing for Granger Causality with Mixed Frequency Data. 
%             Working Paper at UNC. 
% -------------------------------------------------------------------
% Note: Null hypothesis is no causality from the from(1)-th, 
%       (from(1) + iotap)-th, (from(1) + 2*iotap)-th, ... ,
%       from(end)-th variables to to(1)-th, (to(1) + iota)-th, 
%       (to(1) + 2*iota)-th, ... , to(end)-th variables.
%       This code can only deal with constant increments.
%       (e.g., it cannot test no causality from the first, second, 
%              and fourth variables.)
% --------------------------------------------------------------------
% See Also: VAR_est1.m, Wald_test.m, sim_phauto.m
% --------------------------------------------------------------------
% Written by Kaiji Motegi, UNC Chapel Hill.
% Jun. 28, 2013.
% --------------------------------------------------------------------

% error checking
if max(size(from)) > 2
    error('The dimension of "from" must be at most 2.')
end;
if max(size(to)) > 2
    error('The dimension of "to" must be at most 2.')
end;

% Default setting
if nargin == 3
    iota = 1; iotap = 1; bs = 'off';
elseif nargin == 4
    iotap = 1; bs = 'off';
elseif nargin == 5
    bs = 'off';
end;

% The null hypothesis can be written as:
% A_{k}^{(h)}(a:iota:b, c:iotap:d) = zeros(dimto, dimfrom) for k = 1,...,p
% i.e., A_{k}^{(h)}'(c:iotap:d, a:iota:b) = zeros(dimfrom, dimto) for k = 1,...,p
a = to(1);   b = to(end);
c = from(1); d = from(end);

% error checking
if mod(b - a, iota) ~= 0
    error('(b-a) / iota must be a nonnegative integer.');
end;
if mod(d - c, iotap) ~= 0
    error('(d-c) / iotap must be a nonnegative integer.');
end;

dimto = (b - a) / iota + 1;      % dimension of "to" variables
dimfrom = (d - c) / iotap + 1;   % dimension of "from" variables

% load basic quantities
p = result_VAR_est.nlag;   % VAR lag length
q = dimfrom * dimto * p;   % # of zero restrictions
K = result_VAR_est.ndim;   % dimension of VAR process

% construct delta, a key quantity for constructing selection matrix R
% See Eq. (3.5) in GHM
delta = zeros(dimto * p, 1);
delta(1) = (a - 1) * p * K + c;
if dimto*p > 1
    for l = 2:(dimto*p)
         % Basic increment is K.
         delta(l) = delta(l-1) + K;
         % When l-1 is a multiple of p, an extra increment is needed. 
         if mod(l-1,p) == 0
             delta(l) = delta(l) + p * K * (iota - 1);
         end;    
    end;    
end;    

% construct Lambda and selection matrix R
R = zeros(q, p * K^2);
for l = 1:(dimto*p)
     % construct Lambda
     Lambda = zeros(dimfrom, p * K^2);
     for j = 1:dimfrom
          Lambda(j, delta(l)+(j-1)*iotap) = 1;
     end;    
     % construct R
     first = (l-1)*dimfrom + 1;
     last = l * dimfrom;
     R(first:last, :) = Lambda;
end;    

% Compute r
r = zeros(q, 1);
    
% Implement Wald test (without bootstrap)
[W, P_value] = Wald_test(result_VAR_est, R, r);
    
% Bootstrap (optional)
if ischar(bs) == 0   % implement bootstrap
    h = result_VAR_est.horizon;
    nobs = result_VAR_est.T_star;

    if h > 1
         % generate random numbers
         bsrnd = randn(nobs, K, bs);
        
         % import basic quantities
         lambda = result_VAR_est.lambda;
         %Omega1 = result_VAR_est.Omega1;
         %L = chol(Omega1, 'lower');        % Cholesky decomposition
         B1 = result_VAR_est.OLSE1;
         Bh = result_VAR_est.OLSE;
    
         % impose the null hypothesis of non-causality at horizon h
         for i = 1:p
              first = (i-1) * K + 1;
              last = i * K;
              temph = Bh(first:last, :)';
              temph(a:iota:b, c:iotap:d) = 0;
              Bh(first:last, :) = temph';        
         end;
 
         % Compute Wald statistics based on simulated data
         W_bs = zeros(bs, 1);
         for n = 1:bs
              %E = bsrnd(:,:,n) * L';
              E = bsrnd(:, :, n) .* result_VAR_est.residmat;
              X = sim_phauto(E, Bh, B1, p, h);
              result = VAR_est1(X, p, h, lambda);
              % disp('stab cond = '); disp(result.stab);
              % if analytical covariance matrix is used, Sigma needs to be
              % replaced. This is relevant only for benchmark simulation.
              if result_VAR_est.analcov_flag == 1
                  result.Sigma = result_VAR_est.Sigma;
              end;    
              [W_bs(n), ~] = Wald_test(result, R, r);
         end;
         % boostraped p-value
         P_value = (1 + sum(W_bs >= W)) / (bs + 1);    

    elseif h == 1         
        
        %% This is an old part where bootstrap was available for h=1 only.
         lambda = result_VAR_est.lambda;
         stdnorm = randn(nobs, K, bs);
         %Omega1 = result_VAR_est.Omega1;
         %L = chol(Omega1, 'lower');
         B1 = result_VAR_est.OLSE1;
         Amat = zeros(K, K, p);
         for i = 1:p
              first = (i-1) * K + 1;
              last = i * K;
              temp = B1(first:last, :)';
              temp(a:iota:b, c:iotap:d) = 0;
              Amat(:,:,i) = temp;
         end;
         W_bs = zeros(bs, 1);
         for n = 1:bs
              %E = stdnorm(:,:,n) * L';
              E = stdnorm(:,:,n) .* result_VAR_est.residmat;
              X = sim_VAR(E, Amat);
              result = VAR_est1(X, p, h, lambda);
              if result_VAR_est.analcov_flag == 1
                  result.Sigma = result_VAR_est.Sigma;
              end;    
              [W_bs(n), ~] = Wald_test(result, R, r);
         end;
         P_value = (1 + sum(W_bs >= W)) / (bs + 1);    
    end;     
end;    
