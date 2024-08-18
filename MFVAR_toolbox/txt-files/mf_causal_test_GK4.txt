function [W, P_val] = mf_causal_test_GK4(result_VAR_est, m, K_H, casenum, fromind, toind, bs)
% PURPOSE: Implement Granger causality tests based on mixed frequency VAR models
%--------------------------------------------------------------------------
% USAGE: [W, P_val] = mf_causal_test_GK3(result_VAR_est, m, K_H, casenum, fromind, toind, bs)
% where: result_VAR_est = result of VAR_est.m  (structure)
%                     m = # of high frequency observations in one 
%                         unit of low frequency period
%                   K_H = # of high frequency variables
%               casenum = 1 for "low to low"
%                         2 for "high to low"
%                         3 for "low to high"
%                         4 for "high to high"
%               fromind = index for the "from" variable (scalar)
%                         This must be smaller than max{K_H, K_L}
%                         (default for Cases 2 & 3 = 1)
%                 toind = index for the "to" variable (scalar)
%                         This must be smaller than max{K_H, K_L}
%                         (default for Cases 2 & 3 = 1)
%                    bs = # of replications for bootstrap
%                         'off' if bootstrap is not used (default = 'off')
%--------------------------------------------------------------------------
% RETURNS: W     = Wald statistic
%          P_val = Associated P value 
%--------------------------------------------------------------------------
% References: E.Ghysels, J.B.Hill & K.Motegi (2013)
%             Granger Causality in Mixed Frequency Vector 
%             Autoregressive Models. Working Paper at UNC. 
% -------------------------------------------------------------------------
% Note: Null hypothesis for each case is as follows:
%          Case 1. no causality from the fromind-th low-freq variable to the
%                  toind-th low-freq variable.
%          Case 2. no causality from the fromind-th high-freq variable to the
%                  toind-th low-freq variable.
%          Case 3. no causality from the fromind-th low-freq variable to the
%                  toind-th high-freq variable.
%          Case 4. no causality from the fromind-th high-freq variable to the
%                  toind-th high-freq variable.
%       This code can only deal with one "from" variable and one "to" 
%       variable.
%       Also, this code is assuming that all high frequency variables are
%       put before all low frequency variables. Correct order is:
%            [x_{H,1} (tau, 1);
%             x_{H,2} (tau, 1);
%                :
%                :
%             x_{H,1} (tau, m);
%             x_{H,2} (tau, m);
%             x_{L,1} (tau)   ;
%             x_{L,2} (tau)   ];
%       Make sure that the following is a WRONG order:
%            [x_{H,1} (tau, 1);
%                :
%             x_{H,1} (tau, m);
%             x_{H,2} (tau, 1);
%                :
%             x_{H,2} (tau, m);
%             x_{L,1} (tau)   ;
%             x_{L,2} (tau)   ];
% -------------------------------------------------------------------------
% See Also: VAR_est1.m, Wald_test.m, causality_test_GK4.m
% -------------------------------------------------------------------------
% Written by Kaiji Motegi, UNC Chapel Hill.
% Aug. 4, 2012.
% -------------------------------------------------------------------------

% default for Bootstrap flag = 'off'
if nargin < 7
    bs = 'off';
end;    

% dimension of MFVAR process
K = result_VAR_est.ndim;
% # of low frequency variables
K_L = K - m * K_H;

% default increments
iota = 1; iotap = 1;

switch casenum
    case 1   % causality from the fromind-th low-freq variable to the
             % toind-th low-freq variable
             
        % Error checking
        if fromind > K_L
            error('wrong fromind');
        end;
        if toind > K_L
            error('wrong toind');
        end;
        
        a = m * K_H + toind;
        b = a;
        c = m * K_H + fromind;
        d = c;
        
    case 2   % causality from the fromind-th high-freq variable to the
             % toind-th low-freq variable
        
        % Default setting
        if nargin == 4
            fromind = 1;
            toind = 1;
        end;    
             
        % Error checking
        if fromind > K_H
            error('wrong fromind');
        end;
        if toind > K_L
            error('wrong toind');
        end;
             
        a = m * K_H + toind;
        b = a;
        c = fromind;
        d = fromind + (m-1) * K_H;
        iotap = K_H;
        
    case 3   % causality from the fromind-th low-freq variable to the
             % toind-th high-freq variable
             
        % Default setting
        if nargin == 4
            fromind = 1;
            toind = 1;
        end; 
        
        % Error checking
        if fromind > K_L
            error('wrong fromind');
        end;
        if toind > K_H
            error('wrong toind');
        end;     
             
        a = toind;
        b = toind + (m-1) * K_H;
        c = m * K_H + fromind;
        d = c;        
        iota = K_H;
        
    case 4   % causality from the fromind-th high-freq variable to the
             % toind-th high-freq variable
             
        % Error checking
        if fromind > K_H
            error('wrong fromind');
        end;
        if toind > K_H
            error('wrong toind');
        end; 
             
        a = toind;
        b = toind + (m-1) * K_H;
        c = fromind;
        d = fromind + (m-1) * K_H;
        iota = K_H;
        iotap = K_H;
        
    otherwise
        error('Unknown case number');
end;
        
% specify "from" and "to" components
from = [c d];
to = [a b];

% Implement causality test
[W, P_val] = causality_test_GK4(result_VAR_est, from, to, iota, iotap, bs);