%% Weather demo - precipitation and temperature
% Mixed-frequency real-world records:
%   x1, y1: precipitation - monthly (y1 = yearly sum  of x1)
%   x2, y2: temperature   - monthly (y2 = yearly mean of x2)

set(0, 'defaultAxesFontName', 'Calibri')
set(0, 'defaultTextFontName', 'Calibri')

% Figure output folder, relative to the current folder. savefigs needs it
% to be a string and does not create it itself.
fig_folder = "Figures/Weather";
if ~isfolder(fig_folder), mkdir(fig_folder); end

%% Load the precipitation / temperature records

[x1, y1, x2, y2] = precipitation;
x1 = x1'; y1 = y1'; x2 = x2'; y2 = y2';

year_start = 1869; year_end = 2023;       % span covered by precipitation.m
nobs = 1860;                              % 155 years x 12 months
ntrials = 1;
fs = 12;                    % monthly series
fs2 = 1;                    % yearly series
downsample_ratio = 12;

% YOY expects an N x 3 table of (year, month, value).
months = repmat((1:12)', year_end - year_start + 1, 1);
years = repelem((year_start:year_end)', 12);
table_x1 = [years, months, x1(:)];

% %% Percentage change relative to the within-year baseline
% 
% x1_perc = reshape(PercChange(x1, 12)', 1, []);
% figure; plot(x1_perc(1:120))
% 
% x2_perc = reshape(PercChange(x2, 12)', 1, []);
% figure; plot(x2_perc(1:120))

%% Detrend and overview all four series

x1 = detrend(x1); x2 = detrend(x2);
y1 = detrend(y1); y2 = detrend(y2);

figure;
subplot(411); plot(x1, 'LineWidth', 2, 'Color', 'b'); axis tight;
xlabel('Monthly precipitation', 'FontSize', 18)
subplot(412); plot(y1, 'LineWidth', 2, 'Color', 'r'); axis tight;
xlabel('Yearly precipitation', 'FontSize', 18)
subplot(413); plot(x2, 'LineWidth', 2, 'Color', 'b'); axis tight;
xlabel('Monthly temperature', 'FontSize', 18)
subplot(414); plot(y2, 'LineWidth', 2, 'Color', 'r'); axis tight;
xlabel('Yearly temperature', 'FontSize', 18)

%% Monthly vs. yearly precipitation, first 10 years

figure("Position", [0 0 450 300]);
plot(x1(1:120), "LineWidth", 2);
ylabel("Monthly precipitation")
yyaxis right
plot(1:12:121, [nan y1(1:10)], "LineWidth", 3);
ylabel("Yearly precipitation")
ax = gca; set(ax, "fontSize", 20); set(ax, "LineWidth", 3);
axis tight
xticks([1, 61, 121]); xticklabels([0, 5, 10])
xlabel("Year")
% saveas(gcf, fullfile(fig_folder, "Signal_prec"), "epsc")

%% Year-over-year growth of the monthly precipitation
% Computed on the raw series - growth rates of a detrended signal are
% meaningless, since its baseline crosses zero.

yoy_prec = reshape(YOY(table_x1, 12, year_start, year_end)', 1, []);
figure; plot(yoy_prec)

%% Plot 20 years of the yearly precipitation

figure('Position', [100, 100, 400, 250]);
plot(y1(1:20), 'k');
set(gca, "Visible", "off")

%% Granger causality between the monthly series

[GC_f, freq] = G_cause(permute(cat(3, x1, x2), [3,2,1]), fs, ...
    "mode", "tf", "plot", true, "report", true, "maxorder", 30);

%% Spectral GC, precipitation (x1) -> temperature (x2)
% MVGC convention: GC_f(i,j,:) is the causality from variable j to variable i.

figure("Position", [0, 0, 400, 300])
plot(freq, squeeze(GC_f(2,1,:)), "LineWidth", 8);
ax = gca; set(ax, "fontSize", 30); set(ax, "LineWidth", 5);
xlabel(ax, "Freq (cycles/year)")
title(ax, "Prec -> Temp")
ylabel(ax, "GC")
% saveas(gcf, fullfile(fig_folder, "GC_PtoT"), "epsc")

%% Granger causality between the yearly series

[GC_f_yearly, freq_yearly] = G_cause(permute(cat(3, y1, y2), [3,2,1]), fs2, ...
    "mode", "tf", "plot", true, "report", false);

% Same test on the percentage-change series, if the section above that
% builds x1_perc / x2_perc is enabled:
% [GC_f_perc, freq_perc] = G_cause(permute(cat(3, x1_perc, x2_perc), [3,2,1]), fs, ...
%     "mode", "tf", "plot", true, "report", true);

%% Cross-correlation between the monthly precipitation and temperature

[r, l] = xcorr(x1', x2', 100, "normalized");
figure; plot(l, r, "LineWidth", 3)

%% MF-TFCCA parameters
% X is the monthly series (high rate), Y the yearly series (low rate).
% Time is measured in years here, so fs = 12 is "samples per second".

win = 20;                   % STFT window, in month
real_STFT = false;          % false: complex-valued STFT
maxlag1 = 4; maxlag2 = 4;   % lag search range, in years

% The window and the lag range both eat into the usable span, so the
% analysis window is kept clear of both ends of the 155-year record.
time_start = 10; time_end = 145;

surrogate_num = 100;        % single trial, so this drives the runtime

%% MF-TFCCA analysis between monthly prec and yearly temp

tic
[lags, r_0, r_su] = MFTFCCA_lagCC(x1, y2, fs, downsample_ratio, ...
    win, time_start, time_end, maxlag1, maxlag2, real_STFT, ...
    "surrogate_num", surrogate_num, "plot", true, "calibrations", 1, ...
    "partial_lag_set", y2, "n_lags", 5);
toc

savefigs(gca, fig_folder, "high prec to low temp with partial y")

%% MF-TFCCA analysis between yearly prec and monthly temp

tic
[lags_t, r_0_t, r_su_t] = MFTFCCA_lagCC(x2, y1, fs, downsample_ratio, ...
    win, time_start, time_end, maxlag1, maxlag2, real_STFT, ...
    "surrogate_num", surrogate_num, "plot", true, "calibrations", 1, ...
    "partial_lag_set", y1, "n_lags", 5);
toc

savefigs(gca, fig_folder, "high temp to low prec with partial y")


%% MF-TFCCA spectral coefficients at one lag
% Pick the lag from the profile above, e.g. its peak:
%   [~, i_peak] = max(r_0); lag = lags(i_peak);

lag = 0;

tic
MFTFCCA_spectral(x1, y2, fs, downsample_ratio, win, ...
    time_start, time_end, real_STFT, ...
    "lag", lag, "surrogate_num", 20, "plot", true);
toc
