%% Add src folder and all subfolders to path
% addpath(genpath('../src'));

%% Generate bivariate/trivariate VAR system {X1, X2, X3} as X (3*nobs*ntrials)

ntrials   = 100;    % Number of trials
nobs = 1000;    % Number of data points per trial
fs = 200;    % sample rate (Hz)
seed = 2; rng_seed(seed);
dt = 1/fs;
 
%% VAR parameters

AT = zeros(10, 10);
% AT(1:5, 1:5) = 0.15;
% AT(6:10, 6:10) = 0.15;

% TODO: add resonate freq

AT(2,1,1) = -0.4;
AT(2,1,2) = 0.7;
AT(2,1,3) = -0.1;
% 0.63, 0.67, 0.65, 0.60, 0.71, 0.54, 0.63, 0.61, 0.85, 0.69, 0.68, 0.88, 
AT(1,3,10) = 0.63;
AT(1,4,10) = 0.67;
AT(1,5,1) = 0.65;
AT(1,6,1) = 0.60;
AT(1,7,1) = 0.71;
AT(1,8,1) = 0.54;



%% Generate VAR
SIG = diag(ones(10, 1));

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

[S, f, t] = STFT(X(9,:,1)', fs, win, win-downsample_ratio, FFT_num, "fpass", [], "plot", 0);
figure("Position", [0,0,550,300]); imagesc(t, f, mag2db(abs(S(:, 1:fs2*4)))); colormap('jet'); 
cax = colorbar; axis xy; set(cax, 'LineWidth', 3);
set(gca, 'FontSize', 25); set(gca, 'LineWidth', 3);
xticks(t(end)); xticklabels("4s");
title("STFT(X) (db)"); ylabel("Freq (Hz)");

%% Use MF-TFCCA to test the mixed-frequency causality

win = 40;

tic
[lags, r_0, r_su] = MFTFCCA_lagCC(x1(1:end,:), squeeze(permute(X(4,1:downsample_ratio:end,:), [1,3,2])), fs, downsample_ratio, ...
    win, 1, 4, 0.5, 0.5, false, "surrogate_num", 3, "plot", 1 ...
    , "partial_high_set", permute(X([3,5,7],:,:), [1,3,2]), "partial_low_set", squeeze(permute(X([2,6,10],1:downsample_ratio:end,:), [1,3,2])));
toc

% savefigs(gcf, "Figures", "highd-x1-x4-w-partial");
%% Use MFVAR to test the mixed-frequency causality

% X_MF = [reshape(x1(1,:), [nobs/downsample_ratio, downsample_ratio]), x2(1,:)'];
MF_n_trial = 10;
X_MF = cat(1, reshape(x1(1:MF_n_trial,:)', [downsample_ratio, nobs/downsample_ratio, MF_n_trial]), ...
    reshape(x3(1:MF_n_trial,:)', [downsample_ratio, nobs/downsample_ratio, MF_n_trial]), ...
    reshape(x2(1:MF_n_trial,:)', [1, nobs/downsample_ratio, MF_n_trial]));
X_MF = permute(X_MF, [2,1,3]);
size(X_MF)
%%
lag_length = 1;
pval = MFVAR(X_MF, downsample_ratio, 2, 1, lag_length, char('x1','x3','x2'));


%% Use MIDAS to test the mixed-frequency causality
% Syntax:
% [...] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate)
% [...] = MIDAS_ADL(DataY,DataYdate,DataX,DataXdate,name,value,...)

time_stamps_1 = 1:nobs;
time_stamps_2 = 1:downsample_ratio:nobs;
results = MIDAS_ADL(x2(1,:)', time_stamps_2', x1(1,:)', time_stamps_1');
