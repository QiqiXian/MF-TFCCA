
%% Basic Information
% Example code on Ghysels' (2012) mixed frequency vector autoregression.
% We also illustrate how to implement impulse response analysis, forecast error
% variance decomposition, and Ghysels, Hill, and Motegi's (2013)
% mixed frequency Granger causality test.

% Since Ghysels' (2012) MF-VAR share many aspects with classic single-frequency VAR,
% this code should be relevant for those who work on classic VAR as well.

%%%%%%%% Ten required codes in order to run this main code %%%%%%%%%%%%%%%%%%%%%%% 
% 1. VAR_est1.m: Fit (p,h)-autoregression (i.e. VAR(p) model iterated h-times) 
%                with Newey and West's (1987) HAC estimator. Newey and West's (1994) 
%                automatic bandwidth selection is available.
% 2. irf3.m: Compute impulse response function at horizon 0, 1, ..., hmax
%            along with bootstrapped confidence intervals. 
%            Cholesky decomposition is used.
% 3. var_decomp.m: Forecast error variance decomposition at horizon 0, ..., hmax-1.
%                  Cholesky decomposition is used.
% 4. MFCTGK_all1.m: Implement bilateral mixed frequency Granger causality tests
%                   for all possible pairs.
% 5. CTGK_all1.m: Implement bilateral Granger causality tests for all
%                 possible pairs.
% 6. sim_VAR.m: Simulate VAR(p) processes.
% 7. sim_phauto.m: Simulate (p,h)-autoregression (i.e. VAR(p) iterated h-times)
% 8. causality_test_GK4.m: Implement Granger causality tests. 
%                          Goncalves and Killian's (2004) bootstrap is available.
% 9. mf_causal_test_GK4.m: Implement mixed frequency Granger causality tests. 
%                          Goncalves and Killian's (2004) bootstrap is available. 
% 10. Wald_test.m: Implement Wald tests.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Written by Kaiji Motegi, Waseda University
% Sep. 9, 2014.

clear all;
close all;
clc;


%% Scenario
% We analyze three variables x, y, z.
% x is a monthly variable while y and z are quarterly variables.
% The ratio of sampling frequencies, m, is equal to 3.
% construct a 5 x 1 mixed frequency vector X = [x1, x2, x3, y, z]'.
% Assume DGP is MF-VAR(1) with coefficient:
A = [   0,   0.1,   0.4,    0,    0;
      0.1,  -0.1,   0.2,    0,    0;
        0,     0,   0.1,    0,    0;
      0.0,  -0.9,   0.9,  0.2,    0;
        0,     0,     0,  0.9,  0.6];
 
% Evidently, (1) x causes y and (2) y causes z.
% See A(4, 1:3) for (1) and A(5,4) for (2).

% Let's generate data from this DGP and fit MF-VAR(1) to see what happens.


%% Step 1: Initial setting
% general setting
T = 80;            % sample size is 80 quarters, a realistic size.
m = 3;             % ratio of sampling frequencies (month vs. quarter)
K_H = 1;           % one high frequency variable x
K_L = 2;           % two low frequency variables y, z
K = K_L + m*K_H;   % dimension of MF-VAR
p = 1;             % VAR lag length included. 
                   % true lag order is 1.
lambda = 'NW';     % use Newey and West's (1994) automatic bandwidth selection

% Impulse response functions
irfhmax = 6;       % maximum horizon
figureflag = 1;    % draw figure
irfalpha = 0.05;   % draw 95% bootstrapped confidence interval
bsnum = 500;       % # of bootstrap replications
labels = char('x1', 'x2', 'x3', 'y', 'z');  % labels

% forecast error variance decomposition
vdhmax = 6;        % maximum horizon

% Granger causality tests
gcbs = 1999;                        % # of bootstrap replications
dispflag = 1;                       % display p-values
gclabels = char('x', 'y', 'z');     % labels


%% Step 2: Mixed frequency analysis
% generate normal error
E = 0.1 * randn(T, K);
% generate VAR(1) process
Data = sim_VAR(E, A);

% fit MF-VAR(1)
result1 = VAR_est1(Data, p, 1, lambda);
% impulse
[IRF, lb, ub] = irf3(result1, irfhmax, figureflag, irfalpha, bsnum, labels);
% variance decomposition
vd_mf = var_decomp(result1, vdhmax);
% causality test
disp('%%%%% Mixed Frequency, horizon = 1 %%%%%')
pval_mat1 = MFCTGK_all1(result1, m, K_H, gcbs, dispflag, gclabels);
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

result2 = VAR_est1(Data, p, 2, lambda);  
disp('%%%%% Mixed frequency, horizon = 2 %%%%%')
pval_mat2 = MFCTGK_all1(result2, m, K_H, gcbs, dispflag, gclabels);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
disp(blanks(3)');

% ... Now you can see x does cause z at horizon 2.


%% Step 3: Low frequency analysis
% for comparison, aggregate x into quarterly frequency (flow sampling)
Data_ = [mean(Data(:,1:3), 2), Data(:,4), Data(:,5)];
% Data_ = [Data(:,1), Data(:,4), Data(:,5)];

% fit VAR(1) 
result_ = VAR_est1(Data_, p, 1, lambda);
% IRF
[IRF_, lb_, ub_] = irf3(result_, irfhmax, figureflag, irfalpha, bsnum, gclabels);
% variance decomposition
vd_lf = var_decomp(result_, vdhmax);
% causality test
disp('%%%%% Low Frequency, horizon = 1 %%%%%')
pval_mat_lf = CTGK_all1(result_, gcbs, dispflag, gclabels);
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');

%%%%%%%%%%%%%%%%%%%%%% REMARK %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Based on low frequency model, you cannot observe x causing y.
% This is because the positive impact of x3 on y and the negative impact of
% x2 on y offset each other after flow aggregation.
% This highlights an advantage of mixed frequency approach.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
