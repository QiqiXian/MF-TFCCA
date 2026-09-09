%% Generate Logistic Model

nobs = 500;
seed = 1;
ntrials = 100;
fs = 50;
seeds = zeros(ntrials, 1);
x1 = NaN(ntrials, nobs); x2 = NaN(ntrials, nobs);

mode = "l-uni";     % l-uni or l-bi

for i = 1:ntrials
    while sum(~isfinite(x1(i,:))) || sum(~isfinite(x2(i,:)))
        [x1(i, :), x2(i, :)] = generate_nonlinear_system(nobs, seed, mode);
        seed = seed + 1;
    end
    seeds(i) = seed - 1;
end

%% Granger causality test of Logistic system
G_cause(permute(cat(3, x1, x2), [3,2,1]), fs, "mode", "tf", "plot", 1, "report", 1);

%% Down-sample X2
downsample_ratio = 5;
fs2 = fs / downsample_ratio;
x2 = x2(:, 1:downsample_ratio:end);

%% Plot 5s snapshots of X1 (high sampling rate) and X2 (low sampling rate)
time = 5;

figure('Position', [100, 100, 400, 500] ); 
subplot(2, 1, 1); 
plot(x1(1,1:fs*time), 'k');
set(gca, "Visible", "off")
subplot(2, 1, 2);
stairs(x2(1,1:fs2*time), 'k', "LineWidth", 0.8);
set(gca, "Visible", "off")

%% MF-TFCCA test of Logistic system

win = 30;
% add calibration

tic
[lags, r_0, r_su] = MFTFCCA_lagCC(x1, x2, fs, downsample_ratio, ...
    win, 1.5, 8.5, 1, 1, false, "surrogate_num", 3, "plot", 1 ...
    );
toc

%% Use MF-TFCCA to compute spectral on one certain lag

win = 30;
tic
MFTFCCA_spectral(x1(1:end,:), x2(1:end,:), fs, downsample_ratio, win, ...
    1.5, 8.5, false, ...
    "lag", 0.3, "surrogate_num", 3, "plot", false);
toc

%% Use MFVAR to test the mixed-frequency causality

MF_n_trial = 100;
X_MF = cat(1, reshape(x1(1:MF_n_trial,:)', [downsample_ratio, nobs/downsample_ratio, MF_n_trial]), ...
    reshape(x2(1:MF_n_trial,:)', [1, nobs/downsample_ratio, MF_n_trial]));
X_MF = permute(X_MF, [2,1,3]);
size(X_MF)
%%
lag_length = 10;
pval = MFVAR(X_MF, downsample_ratio, 1, 1, lag_length, char('x', 'y'));


%% Generate coupled RL system
nobs = 1000;
seed = 1;
ntrials = 100;
fs = 100;

fig_folder = ".\Figures\0-Paper_figures\Nonlinear-systems\RL\";

downsample_ratio = 5;
fs2 = fs / downsample_ratio;

seeds = zeros(ntrials, 1);
X = NaN(3, ntrials, nobs); Y = NaN(3, ntrials, nobs);

for i = 1:ntrials
    while sum(~isfinite(X(:,i,:)), 'all') || sum(~isfinite(Y(:,i,:)), 'all')
        [X(:, i, :), Y(:, i, :)] = generate_nonlinear_system(nobs, seed, "RL");
        seed = seed + 1;
    end
    seeds(i) = seed - 1;
    
end

% X = X + rand(size(X)) *.5; Y = Y + rand(size(Y)) *.5;

x1 = squeeze(X(1,:,:)); x2 = squeeze(X(2,:,:)); x3 = squeeze(X(3,:,:));
y1 = squeeze(Y(1,:,:)); y2 = squeeze(Y(2,:,:)); y3 = squeeze(Y(3,:,:));

%% GC test on RL system
[f,~] = G_cause(permute(cat(1,X,Y), [1,3,2]), fs, "mode", 2, "plot", 1, "report", 0, "seed", 5);

%% Downsample

%% GC test on the high-sampling-rate part
G_cause(permute(X, [1,3,2]), fs, "mode", 2, "plot", 1, "report", 0);
