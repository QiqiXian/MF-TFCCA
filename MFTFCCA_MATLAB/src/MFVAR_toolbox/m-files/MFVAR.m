function pval_mat1 = MFVAR(Data, m, K_H, K_L, p, gclabels)
% Data: Time * K * #trials
% Data, m(ratio of sampling frequencies), K_H, K_L, p(VAR lag length
% included), labels, gclabels

%% Step 1: Initial setting
% general setting
% T = 80;            % sample size is 80 quarters, a realistic size.
% m = 3;             % ratio of sampling frequencies (month vs. quarter)
% K_H = 1;           % one high frequency variable x
% K_L = 2;           % two low frequency variables y, z
K = K_L + m*K_H;   % dimension of MF-VAR
% p = 1;             % VAR lag length included. 
%                    % true lag order is 1.
lambda = 'NW';     % use Newey and West's (1994) automatic bandwidth selection

% Impulse response functions
irfhmax = 6;       % maximum horizon
figureflag = 1;    % draw figure
irfalpha = 0.05;   % draw 95% bootstrapped confidence interval
bsnum = 500;       % # of bootstrap replications
% labels = char('x1', 'x2', 'x3', 'y', 'z');  % labels

% forecast error variance decomposition
vdhmax = 6;        % maximum horizon

% Granger causality tests
gcbs = 100;                        % # of bootstrap replications % 1999
dispflag = 1;                       % display p-values
% gclabels = char('x', 'y', 'z');     % labels


%% Step 2: Mixed frequency analysis
% generate normal error
% E = 0.1 * randn(T, K);
% generate VAR(1) process
% Data = sim_VAR(E, A);

% fit MF-VAR(1)
tic
result1 = VAR_est1(Data, p, 1, lambda);
toc
% impulse
% [IRF, lb, ub] = irf3(result1, irfhmax, figureflag, irfalpha, bsnum, labels);
% variance decomposition
% vd_mf = var_decomp(result1, vdhmax);
% causality test
disp('%%%%% Mixed Frequency, horizon = 1 %%%%%')
tic
pval_mat1 = MFCTGK_all1(result1, m, K_H, gcbs, dispflag, gclabels);
toc
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(blanks(3)');

%%%%%%%%%%%%%%%% REMARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ... What can you tell from these results?
% Causality test says (1) x causes y, (2) y causes z, and there is no other causality.
%
% IRF verifies (1) and (2). IRF gives you one more important implication, though. 
% It seems that x does have a significant impact on z! How is this ever possible?
%
% ... This is a typical example of "causal chain". x does cause z via y.
% To see this point, note that:
%
% A^2 = [ 0.01,  -0.01,  0.06,    0,    0;
%        -0.01,   0.02,  0.04,    0,    0;
%            0,      0,  0.01,    0,    0; 
%        -0.09,  -0.09,  0.09, 0.04,    0;
%            0,  -0.81,  0.81, 0.72, 0.36]
%
% The lower-left block is no longer zeros.

% In bivariate case causal chains are never possible, but in more general
% cases causal chains are of great importance.
% To capture causality from x to z, we need to run two-step-ahead causality test.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% result2 = VAR_est1(Data, p, 2, lambda);  
% disp('%%%%% Mixed frequency, horizon = 2 %%%%%')
% pval_mat2 = MFCTGK_all1(result2, m, K_H, gcbs, dispflag, gclabels);
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
% disp(blanks(3)');

% ... Now you can see x does cause z at horizon 2.


% %% Step 3: Low frequency analysis
% % for comparison, aggregate x into quarterly frequency (flow sampling)
% Data_ = [mean(Data(:,1:K_H+1), 2), Data(:,K_H+2)];
% 
% % fit VAR(1) 
% result_ = VAR_est1(Data_, p, 1, lambda);
% % IRF
% [IRF_, lb_, ub_] = irf3(result_, irfhmax, figureflag, irfalpha, bsnum, gclabels);
% % variance decomposition
% vd_lf = var_decomp(result_, vdhmax);
% % causality test
% disp('%%%%% Low Frequency, horizon = 1 %%%%%')
% pval_mat_lf = CTGK_all1(result_, gcbs, dispflag, gclabels);
% disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%%%%%%%%%%%%%%%%%%%%%% REMARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on low frequency model, you cannot observe x causing y.
% This is because the positive impact of x3 on y and the negative impact of
% x2 on y offset each other after flow aggregation.
% This highlights an advantage of mixed frequency approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
