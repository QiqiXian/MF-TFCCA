function [W, P_value] = Wald_test(result_VAR_est, R, r)
% PURPOSE: Implement Wald test for VAR models
%--------------------------------------------------------------
% USAGE: [W, P_value] = Wald_test(result_VAR_est, R, r)
% where: result_VAR_est = result of VAR_est.m (structure)
%                     R = restriction (q-by-pK^2)
%                     r = restriction (q-by-1)
%                     q = # of restrictions
%                     p = # of lags
%                     K = dimension of VAR
%--------------------------------------------------------------
% RETURNS: W       = Wald statistic
%          P_value = Associated P value 
%--------------------------------------------------------------
% References: E.Ghysels, J.B.Hill & K.Motegi (2013)
%             Granger Causality in Mixed Frequency Vector 
%             Autoregressive Models. Working Paper at UNC. 
% --------------------------------------------------------------
% Note: Null hypothesis is R * vec(B(h)) = r, where OLSE for B(h)
%       is saved in result_VAR_est.OLSE.
% ---------------------------------------------------------------
% See Also: VAR_est1.m
% ---------------------------------------------------------------
% Written by Kaiji Motegi, UNC Chapel Hill.
% Aug. 4, 2012.
% ---------------------------------------------------------------

% # of restrictions
q = size(R, 1);

% Error checking
if size(R, 2) ~= result_VAR_est.nlag * (result_VAR_est.ndim)^2
    error('Wrong number of columns for R');
end;
if size(r, 1) ~= q
    error('Wrong number of rows for r');
end;

% Wald test statistic
% See Eq. (2.9) in GHM
%temp1 = R * vec(result_VAR_est.OLSE) - r;
p = result_VAR_est.nlag;
K = result_VAR_est.ndim;
temp1 = R * reshape(result_VAR_est.OLSE, p*K^2, 1) - r;
temp2 = R * result_VAR_est.Sigma * R';
W = result_VAR_est.nobs * temp1' * (temp2 \ temp1);

% Associated P value
% Asymptotic distribution under null hypothesis is chi-squared(q)
% See Theorem 2.2
P_value = 1 - chi2cdf(W, q);
