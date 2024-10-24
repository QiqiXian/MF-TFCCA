%% Generate bivariate/trivariate VAR system {X1, X2, X3} as X (3*nobs*ntrials)

ntrials   = 100;    % Number of trials
nobs = 4000;    % Number of data points per trial
fs = 200;    % sample rate (Hz)
seed = 2; rng_seed(seed);
dt = 1/fs;

%% Resonating frequencies parameters

% f1, f2 for X1; f3, f4 for X2; f5 for X3. 
% Set the corresponding r as 0 to ignore the frequency.

f1 = 80; f2 = 4; f3 = 15; f4 = 2; f5 = 30; 
r1 = 0.9; r2=0.8; r3 = 0.85; r4 = 0.7; r5 = 0.9;

% f1 = 75; f2 = 4; f3 = 15; f4 = 2; f5 = 30; 
% r1 = 0.; r2=0.8; r3 = 0.9; r4 = 0.; r5 = 0.9;

% f1 = 80; f2 = 5; f3 = 10; f4 = 2; f5 = 15; 
% r1 = 0.95; r2=0.7; r3 = 0.85; r4 = 0.7; r5 = 0.8;

%%
% Generate the diagonal elements for VAR coefficient matrix AT
theta1 = f1 * dt * 2 * pi; cos1 = cos(theta1);
theta12 = f2 * dt * 2 * pi; cos12 = cos(theta12);
theta2 = f3 * dt * 2 * pi; cos2 = cos(theta2);
theta22 = f4 * dt * 2 * pi; cos22 = cos(theta22); 
theta3 = f5 * dt * 2 * pi;

clear AT

AT(:,:,1) = [2*(r1*cos1+r2*cos12), 0.0, 0;
             0, 2*r3*cos2 + 2*r4*cos22, 0;
             -0, 0, 2*r5*cos(theta3)];
AT(:,:,2) = [-r1^2-r2^2-4*r1*r2*cos1*cos12, 0, 0;
             0, -r3^2-r4^2-4*r3*r4*cos2*cos22, 0;
             0., 0, -r5^2];
AT(:,:,3) = [2*r1*r2*(r1*cos12+r2*cos1) 0.0 0;
             0.0 2*r3*r4*(r3*cos2+r4*cos22) 0;   % 0.01, 0.5
             -0 0 0];
AT(:,:,4) = [-r1^2*r2^2 -0.0 0;
             0. -r3^2*r4^2 0;
             0 0 0];
         
%% Causality parameters
% Set the non-diagonal elements of AT to control the causality
% AT(i,j,lag) controls the causality from Xj to Xi at lag

% System1: X to Y
AT(2,1,1) = -0.4;
AT(2,1,2) = 0.7;
AT(2,1,3) = -0.1;

% % System2: Y to X
% AT(1,2,1) = 0.05;
% AT(1,2,2) = -0.05;
% AT(1,2,3) = 0.1;

% % System3: Bidirectional
% b = [-0.35, 0.7, -0.3] / 2;
% for i = 1:length(b)
%     AT(1,2,i+20) = b(i);
% end
% b = fir1(10, [0.2, 0.6], hann(11)) / 5;
% for i = 1:length(b)
%     AT(2,1,i+31) = b(i);
% end

% % System4: Chain
% % X1 -> X2
% AT(3,1,1) = -0.7;
% AT(3,1,2) = 1.5;
% AT(3,1,3) = 1;
% % X2 -> Y
% AT(2,3,1) = -0.3;
% AT(2,3,2) = 0.4;
% AT(2,3,3) = -0.3;

% % System5: Parallel
% % X1 -> Y
% AT(2,1,1) = -0.4;
% AT(2,1,2) = 0.7;
% AT(2,1,3) = -0.1;
% % X2 -> Y
% AT(2,3,1) = -0.3;
% AT(2,3,2) = 0.4;
% AT(2,3,3) = -0.3;

%% Generate VAR
SIG = [1 0 0; 0 1 0; 0 0 1];
X = var_to_tsdata(AT, SIG, nobs, ntrials);

x1 = squeeze(X(1,:,:))';
x2 = squeeze(X(2,:,:))';
x3 = squeeze(X(3,:,:))';


%% Granger causality test on the original system
G_cause(X(:,:,:), fs, "mode", "tf", "plot", 1, "report", 1, "maxorder", 10);

%% Down-sample X2
downsample_ratio = 5;
fs2 = fs / downsample_ratio;
x2 = x2(:, 1:downsample_ratio:end);

%% Create phase-amplitude coupled (PAC) X1

f_a = 90;
x1 = x1 - min(x1, [], 2);
x1 = x1 .* (sin(2*pi*f_a*(1:nobs)/fs));

%% Plot 2s snapshots of X1 (high sampling rate) and X2 (low sampling rate)
figure('Position', [100, 100, 400, 500] ); 
subplot(2, 1, 1); 
plot(x1(1,1:fs*2), 'k');
set(gca, "Visible", "off")
subplot(2, 1, 2);
plot(x2(1,1:fs2*2), 'k');
set(gca, "Visible", "off")

figure; plot(x3(1,1:fs*2), 'k')

%% Compute and plot of X1 spectrum sample

win = 30; % window length, number of X1 data points within a window
FFT_num = win; % FFT number

[S, f, t] = STFT(x1(1,:)', fs, win, win-downsample_ratio, FFT_num, "fpass", [], "plot", 0);
figure("Position", [0,0,550,300]); imagesc(t, f, mag2db(abs(S(:, 1:fs2*4)))); colormap('jet'); 
cax = colorbar; axis xy; set(cax, 'LineWidth', 3);
set(gca, 'FontSize', 25); set(gca, 'LineWidth', 3);
xticks(t(end)); xticklabels("4s");
title("STFT(X) (db)"); ylabel("Freq (Hz)");

%% Use MF-TFCCA to test the mixed-frequency causality

win = 40;

tic
[lags, r_0, r_su] = MFTFCCA_lagCC(x1(1:10,:), x2(1:10,:), fs, downsample_ratio, ...
    win, 2, 18, 0.5, 0.5, false, [], "surrogate_num", 300, "plot", 0);
toc

%%
win = 40;
tic
MFTFCCA_spectral(x1(2,:), x2(2,:), fs, downsample_ratio, win, 0.08, 19.9, false, "lag", 0, "surrogate_num", 100, "plot", false);
toc

%% Use MFVAR to test the mixed-frequency causality

% X_MF = [reshape(x1(1,:), [nobs/downsample_ratio, downsample_ratio]), x2(1,:)'];
MF_n_trial = 1;
X_MF = cat(1, reshape(x1(1:MF_n_trial,:)', [downsample_ratio, nobs/downsample_ratio, MF_n_trial]), ...
    reshape(x2(1:MF_n_trial,:)', [1, nobs/downsample_ratio, MF_n_trial]));
    %     reshape(x3(1:MF_n_trial,:)', [downsample_ratio, nobs/downsample_ratio, MF_n_trial]), ...
X_MF = permute(X_MF, [2,1,3]);
size(X_MF)
%%
lag_length = 1;
pval = MFVAR(X_MF, downsample_ratio, 1, 1, lag_length, char('x1','x2','x3','x4','x5',...
    'y'), char('x','y'));