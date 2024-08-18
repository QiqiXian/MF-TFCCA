function pval_mat = MFCTGK_all1(result_VAR_est, m, K_H, bs, dispflag, labels)
% PURPOSE: Implement MFVAR-based Granger causality test with Goncalves and
%          Killian's (2004) parametric bootstrap for all pairs of variables
%---------------------------------------------------------------------
% USAGE: pval_mat = MFCTGK_all(result_VAR_est, m, K_H, bs, dispflag, labels)
% where: result_VAR_est = result of VAR_est1.m  (structure)
%        m   = Ratio of sampling frequencies (scalar)
%        K_H = # of high frequency variables (scalar)
%        bs  = # of bootstrap replications  (scalar)
%              'off' if you do not use bootstrap
%              (default = 'off')
%        dispflag = 1 if you display p-values for each pair
%        labels   = K_star x 1 strings cotaining the name of each variable
%                   e.g. labels = ['SP ';
%                                  'I  ';
%                                  'FFR'];
%                   (default = ['hfvar1'  ;
%                               'hfvar2'  ;
%                                ...
%                               'hfvarK_H';
%                               'lfvar1'  ;
%                               'lfvar2'  ;
%                                ...
%                               'lfvarK_L']; )
%        K_star = K_H + K_L
%        K_L    = # of low frequency variables
%----------------------------------------------------------------------
% RETURNS: pval_mat = matrix containing p-values for all pairs (K_star x K_star)
%                     pval_mat(i,j) has a p-value for the null hypothesis
%                     that variable j does not cause variable i.
%                     Note that variables are ordered as [hfvar1, ..., hfvarK_H,
%                     lfvar1, ..., lfvarK_L].
% ----------------------------------------------------------------------
% Note: This code is assuming that all high frequency variables are
%       put before all low frequency variables when MF-VAR is fitted. 
%       Correct order for the mixed frequency vector X is:
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
% ------------------------------------------------------------------------
% References: E.Ghysels, J.B.Hill, & K.Motegi (2013)
%             Testing for Granger Causality with Mixed Frequency Data. 
%             Working Paper at UNC. 
% ---------------------------------------------------------------------
% See Also: VAR_est1.m, mf_causality_test_GK4.m
% --------------------------------------------------------------------
% Written by Kaiji Motegi, Waseda University.
% Aug. 18, 2014.
% ----------------------------------------------------------------------

K = result_VAR_est.ndim;  % dimension of MF-VAR
K_L = K - m * K_H;        % # of low frequency variables
K_star = K_H + K_L;       % # of variables 

% default setting
if nargin == 3
    bs = 'off';     % no bootstrap
    dispflag = 0;   % do not display results
elseif nargin == 4
    dispflag = 0;
end;

if (dispflag == 1) && (nargin == 5)   % use default labels
    % high frequency variables
    index_HF = (1:K_H)';
    letters_HF = repmat('hfvar', K_H, 1);
    labels_HF = strcat(letters_HF, num2str(index_HF));   % combine 'hfvar' and number
    
    % low frequency variables
    index_LF = (1:K_L)';
    letters_LF = repmat('lfvar', K_L, 1);
    labels_LF = strcat(letters_LF, num2str(index_LF));   % combine 'lfvar' and number 
    
    % combine high frequency variables and low frequency variables
    labels = char(labels_HF, labels_LF);
end;    
         
% Below we consider non-causality from variable j to variable i.
% When i and j are fixed, casenum is determined as follows.
casenum_UL = 4 * ones(K_H, K_H);   % high-to-high causality
casenum_UR = 3 * ones(K_H, K_L);   % low-to-high causality
casenum_LL = 2 * ones(K_L, K_H);   % high-to-low causality
casenum_LR = ones(K_L, K_L);       % low-to-low causality
casenum_mat = [casenum_UL, casenum_UR;
               casenum_LL, casenum_LR];

% fromind is determined as follows.
fromind_mat = [repmat(1:K_H, K_star, 1), repmat(1:K_L, K_star, 1)];
% toind is determined as follows.
toind_mat = [repmat((1:K_H)', 1, K_star); repmat((1:K_L)', 1, K_star)];


% consider all possible pairs of variables
pval_mat = zeros(K_star, K_star);
for i = 1:K_star
     for j = 1:K_star    % test for non-causality from variable j to variable i
          if j ~= i    % do nothing if j = i
              % pick proper casenum, fromind, and toind
              [~, pval_mat(i,j)] = mf_causal_test_GK4(result_VAR_est, m, ...
                                        K_H, casenum_mat(i,j), fromind_mat(i,j), toind_mat(i,j), bs);
              
              if dispflag == 1  % display p-value if desired
                  % strcat deletes blanks while [ ] does not.
                  sentence = ['H_0: ', strcat(labels(j,:)), ' does not cause ', strcat(labels(i,:))];
                  disp(sentence);
                  pval = ['   p-value = ', num2str(pval_mat(i,j))];
                  disp(pval);
                  disp(blanks(1));  % single blank line
              end;    
          end;    
     end;    
end;    
