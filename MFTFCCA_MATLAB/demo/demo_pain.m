%% Pain demo - ACC and S1 recordings
% Laser-evoked recordings, 2 x 2000 x ntrials at fs = 200 Hz.
%   x1: ACC (high sampling rate)
%   x2: S1  (down-sampled to form the low-rate series)
% Only the second half of each trial (samples 1001:end) is analysed.

set(0, 'defaultAxesFontName', 'Calibri')
set(0, 'defaultTextFontName', 'Calibri')

% Figure output folder, relative to the current folder.
fig_folder = "Figures/Pain";
if ~isfolder(fig_folder), mkdir(fig_folder); end

%% Load the recordings

load(fullfile("..", "..", "Data", "ACC_S1.mat"), "data");

fs = 200;
downsample_ratio = 2;
ntrials = size(data, 3);
nobs = 1000;

x1 = squeeze(data(1, 1001:end, :))';     % ACC
x2 = squeeze(data(2, 1001:end, :))';     % S1

%% Granger causality, per trial and pooled

f_all_1 = []; f_all_2 = [];
for i = 1:ntrials
    [f, ff] = G_cause(data(:, 1001:end, i), fs, ...
        "mode", "tf", "plot", false, "report", false, "seed", 1, "maxorder", 50);
    f_all_1 = [f_all_1; f(1,2,:)];
    f_all_2 = [f_all_2; f(2,1,:)];
end

[f, ff] = G_cause(data(:, 1001:end, :), fs, ...
    "mode", "tf", "plot", true, "report", true, "seed", 1, "maxorder", 15);

%% Spectral GC

figure("Position", [0,0,400,300])
plot(ff, squeeze(f(1,2,:)), "LineWidth", 8);
ax = gca; set(ax, "fontSize", 30); set(ax, "LineWidth", 5);
xlabel(ax, "Freq (Hz)")
ylabel(ax, "GC")
axis tight
title("ACC -> S1")
% saveas(gcf, fullfile(fig_folder, "GC_AtoS"), "png")

%% Plot 1s of S1, original vs. down-sampled

figure('Position', [100, 100, 400, 250] );
plot(x2(1,1:1:fs), 'k'); hold on; set(gca, "LineWidth", 2);
stairs(1:downsample_ratio:fs, x2(1,1:downsample_ratio:fs), 'b');
set(gca, 'FontSize', 18);
set(gca, "Visible", "off")
% legend(["Original", "Downsampled"]);
saveas(gcf, fullfile(fig_folder, "signal_S1_downsampled"), "png")

%% Mixed-frequency cross-correlation, both directions

MF_xcorr(x1, x2(:,1:downsample_ratio:end), fs, downsample_ratio, 0.1, 0.1)
saveas(gcf, fullfile(fig_folder, "XCorr_high_ACC_low_S1"), "png")

MF_xcorr(x2, x1(:,1:downsample_ratio:end), fs, downsample_ratio, 0.1, 0.1)
saveas(gcf, fullfile(fig_folder, "XCorr_low_ACC_high_S1"), "png")

%% MF-TFCCA parameters and calibrations

win = 40;
real_STFT = true;
time_start = 0.8; time_end = 4.2;       % seconds
maxlag1 = 0.5; maxlag2 = 0.5;           % seconds

cals = MFTFCCA_compute_calibrations(x2, fs, downsample_ratio, ...
    win, time_start, time_end, maxlag1, maxlag2, real_STFT, "plot", true);
figure; plot(cals(:,1), '.'); hold on; plot(cals(:,2), '.')

%% MF-TFCCA lag-CC profile, high-rate S1 -> low-rate ACC

[lags, r_0, r_su] = MFTFCCA_lagCC(x2, x1(:,1:downsample_ratio:end), fs, downsample_ratio, ...
    win, time_start, time_end, maxlag1, maxlag2, real_STFT, ...
    "calibrations", cals, "surrogate_num", 3);
% saveas(gcf, fullfile(fig_folder, "high_S1_low_ACC-lag-corr-win40"), "png")

%% MF-TFCCA spectral coefficients, high-rate ACC -> low-rate S1

[stft_f, r_0, As0, r_su, ~] = MFTFCCA_spectral(x1, x2(:,1:downsample_ratio:end), fs, downsample_ratio, ...
    win, time_start, time_end, false, "calibrations", cals, "lag", -0.13, ...
    "plot", true, "plot_interval", 0, "surrogate_num", 3, "fpass", []);
set(gca, "YLim", [0 3]);
mean(r_0)
std(r_0)
% saveas(gcf, fullfile(fig_folder, "high_ACC_low_S1-spectral-AtoS"), "png")

%% Band-stop sweep - filter one band out of S1 at a time

bands = [[0.1, 20]; [10,30]; [20, 40]; [30,50]; [40,60]; [50, 70]; [60,80]]';
r_fil = []; r_fil_su = [];
for band = bands
    x2_filtered = x2;
    for i = 1:ntrials
        x2_filtered(i,:) = bandstop(x2(i, :), band, fs);
    end

    [~, r_fil1, ~, r_fil1_su, ~] = MFTFCCA_spectral(x2_filtered, x1(:,1:downsample_ratio:end), ...
        fs, downsample_ratio, win, time_start, time_end, false, ...
        "calibrations", [], "lag", 0.1075, "plot", true, "plot_interval", 0, ...
        "surrogate_num", 5, "fpass", []);
    set(gca, "YLim", [0 3]);

    r_fil = [r_fil; r_fil1];
    r_fil_su = [r_fil_su; mean(r_fil1_su,2)'];
end

%% CC gain relative to the unfiltered signal

figure("Position",  [0 0 400 500])
tiledlayout(5,1, 'TileSpacing', "compact")

nexttile([3,1])
ax = gca; hold on;
set(ax, "fontSize", 25); set(ax, "LineWidth", 4)
errorbar(ax, -1 , 0, std(r_0) / sqrt(ntrials), '.', "MarkerSize", 15, "LineWidth", 4)
errorbar(ax, 0:1:size(r_fil,1)-1, mean(r_fil, 2) - mean(r_0), ...
    std(r_fil, [], 2) / sqrt(ntrials), '.', "MarkerSize", 15, "LineWidth", 4)
plot(ax, xlim, [0 0], "--", "Color", "#AAAAAA", "LineWidth", 2)
box off
set(ax, "TickDir", "out")
hold(ax, "off")
ax.XAxis.Visible = 'off';
title("CC gain")
ylabel("Original")
axis padded
ylim([-0.18 0.06])

nexttile([2,1])
ax = gca;
set(ax, "fontSize", 25); set(ax, "LineWidth", 4); hold(ax, "on")
errorbar(ax, -1, 0 , std(r_su(:)) / sqrt(ntrials), '.', "MarkerSize", 15, "LineWidth", 4)
errorbar(ax, 0:1:size(r_fil,1)-1, mean(r_fil_su, 2) - mean(r_su(:)), ...
    std(r_fil_su,[],2) / sqrt(ntrials), '.', "MarkerSize", 15, "LineWidth", 4)
plot(ax, xlim, [0 0], "--", "Color", "#AAAAAA", "LineWidth", 2)
yticks([min(yticks), 0])
ylabel("Surrogate")
xlabel("Filtered-out band")
box off
set(ax, "TickDir", "out")
axis padded
xticks([-1, 0:size(r_fil,1)-1]);
xticklabels(["Original", "1-20Hz", "10-30Hz", "20-40Hz", "30-50Hz", "40-60Hz", "50-70Hz", "60-80Hz"]);
xtickangle(45)
% saveas(gcf, fullfile(fig_folder, "2fold_high_S1_low_ACC_CCgain_AtoS"), "png")

%% MF-TFCCA spectral at a shorter window, high-rate ACC -> low-rate S1

win_short = 20;

[stft_f, ~, As0, ~, ~] = MFTFCCA_spectral(x1, x2(:,1:downsample_ratio:end), fs, downsample_ratio, ...
    win_short, time_start, time_end, false, "calibrations", cals, "lag", -0.07, ...
    "plot", true, "plot_interval", 0, "surrogate_num", 3);
set(gca, "YLim", [0 3]);
saveas(gcf, fullfile(fig_folder, "high_ACC_low_S1-spectral-AtoS"), "png")

%% Reverse direction - low-rate ACC, high-rate S1

cals = MFTFCCA_compute_calibrations(x2, fs, downsample_ratio, ...
    win_short, time_start, time_end, maxlag1, maxlag2, false, "plot", true);
figure; plot(cals(:,1), '.'); hold on; plot(cals(:,2), '.')

[lags, r_0, r_su] = MFTFCCA_lagCC(x2, x1(:,1:downsample_ratio:end), fs, downsample_ratio, ...
    win_short, time_start, time_end, maxlag1, maxlag2, false, ...
    "calibrations", cals, "surrogate_num", 3);
saveas(gcf, fullfile(fig_folder, "high_S1_low_ACC-lag-corr"), "png")

[stft_f, ~, As1, ~, ~] = MFTFCCA_spectral(x2, x1(:,1:downsample_ratio:end), fs, downsample_ratio, ...
    win_short, time_start, time_end, false, "calibrations", [], "lag", -0.03, ...
    "plot", true, "plot_interval", 0, "surrogate_num", 3);
set(gca, "YLim", [0 3]);

% Compare the coefficients at two lags:
% [stft_f, ~, As2, ~, ~] = MFTFCCA_spectral(x2, x1(:,1:downsample_ratio:end), fs, downsample_ratio, ...
%     win_short, time_start, time_end, false, "calibrations", [], "lag", 0.07, ...
%     "plot", true, "plot_interval", 0, "surrogate_num", 3);
% plot_CCA_coeff(stft_f, As1, "As2", As2)
% set(gca, "YLim", [0 3])
% saveas(gcf, fullfile(fig_folder, "high_S1_low_ACC-spectral-StoA"), "png")

[stft_f, ~, As0, ~, ~] = MFTFCCA_spectral(x2, x1(:,1:downsample_ratio:end), fs, downsample_ratio, ...
    win_short, time_start, time_end, false, "calibrations", [], "lag", 0.07, ...
    "plot", true, "plot_interval", 0, "surrogate_num", 3);
set(gca, "YLim", [0 3]);
saveas(gcf, fullfile(fig_folder, "high_S1_low_ACC-spectral-AtoS"), "png")

%% MF-TFCCA runtime vs. number of trials

for MF_n_trial = [1, 10, ntrials]
    tic
    MFTFCCA_spectral(x2(1:MF_n_trial,:), x1(1:MF_n_trial,1:downsample_ratio:end), ...
        fs, downsample_ratio, win_short, time_start, time_end, false, ...
        "calibrations", cals, "lag", 0.10, "plot", false, "plot_interval", 1, ...
        "surrogate_num", 100);
    toc
end

%% MFVAR runtime vs. number of trials
% Recorded timings, 100 surrogates:
%   x2, x1; p=5; trials 1-10:  15.99s; 20-29: 169.38s; 30-39: 145.83s
%   x2, x1; p=1; trials 1-10:   6.17s; 30-39:   7.14s
%   x1, x2; p=1; trials 1-10:  16.76s; 10-19:  17.19s
%   x1, x2; p=5; trials 1-10: 258.00s; 11-20: 210.99s

lag_length = 5;

for MF_n_trial = [1, 10, ntrials]
    X_MF = cat(1, reshape(x1(1:MF_n_trial,:)', [downsample_ratio, nobs/downsample_ratio, MF_n_trial]), ...
        reshape(x2(1:MF_n_trial,1:downsample_ratio:end)', [1, nobs/downsample_ratio, MF_n_trial]));
    X_MF = permute(X_MF, [2,1,3]);

    tic
    pval = MFVAR(X_MF, downsample_ratio, 1, 1, lag_length, char('x','y'));
    toc
end
