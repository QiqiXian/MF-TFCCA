function pval_mat = CTGK_all1(result_VAR_est, bs, dispflag, labels)
% PURPOSE: Implement VAR-based Granger causality test with Goncalves and
%          Killian's (2004) parametric bootstrap for all pairs of variables
%---------------------------------------------------------------------
% USAGE: pval_mat = CTGK_all(result_VAR_est, bs, dispflag, labels)
% where: result_VAR_est = result of VAR_est.m  (structure)
%        bs             = # of bootstrap replications
%                         'off' if you do not use bootstrap
%                         (default = 'off')
%        dispflag = 1 if you display p-values for each pair
%        labels   = K x 1 strings cotaining the name of each variable
%                   e.g. labels = ['SP ';
%                                  'I  ';
%                                  'FFR'];
%                   (default = ['var1';
%                               'var2';
%                                ...
%                               'varK']; )
%        K = dimension of VAR
%----------------------------------------------------------------------
% RETURNS: pval_mat = matrix containing p-values for all pairs (K x K)
%                     pval_mat(i,j) has a p-value for the null hypothesis
%                     that variable j does not cause variable i.
% ----------------------------------------------------------------------
% References: E.Ghysels, J.B.Hill, & K.Motegi (2013)
%             Testing for Granger Causality with Mixed Frequency Data. 
%             Working Paper at UNC. 
% ---------------------------------------------------------------------
% See Also: VAR_est1.m, causality_test_GK4.m
% --------------------------------------------------------------------
% Written by Kaiji Motegi, Waseda University.
% Aug. 18, 2014.
% ----------------------------------------------------------------------

% dimension of VAR
K = result_VAR_est.ndim;

% default setting
if nargin == 1
    bs = 'off';     % no bootstrap
    dispflag = 0;   % do not display results
elseif nargin == 2
    dispflag = 0;
end;

if (dispflag == 1) && (nargin == 3)   % use default labels
    index = (1:K)';
    letters = repmat('var', K, 1);
    labels = strcat(letters, num2str(index));   % combine 'var' and number
end;    
         
% consider all possible pairs of variables
pval_mat = zeros(K, K);
for i = 1:K
     for j = 1:K    % test non-causality from variable j to variable i
          if j ~= i    % do nothing if j = i
              [~, pval_mat(i,j)] = causality_test_GK4(result_VAR_est, j, i, 1, 1, bs);
              
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
