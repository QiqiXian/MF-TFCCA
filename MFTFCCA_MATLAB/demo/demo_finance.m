%% Finance demo
% Mixed-frequency macroeconomic series, year-over-year growth rates:
%   CPI  - monthly   (fs  = 12)
%   OIL  - monthly   (fs  = 12)
%   GDP  - quarterly (fs2 = 4)
% Common span: 1987-2016.

set(0, 'defaultAxesFontName', 'Calibri')
set(0, 'defaultTextFontName', 'Calibri')

% Locate this file, so the demo runs whatever the current folder is.
% "Run Section" evaluates through a temporary copy of the script, so
% mfilename can report a path under tempdir. Try the Editor API and pwd as
% well, and keep only a candidate that really has the data folder below it.
candidates = {};
this_file = mfilename('fullpath');
if ~isempty(this_file), candidates{end+1} = fileparts(this_file); end
try
    candidates{end+1} = fileparts(matlab.desktop.editor.getActiveFilename);
catch
end
candidates{end+1} = pwd;

demo_dir = '';
for i = 1:numel(candidates)
    if ~isempty(candidates{i}) && ...
            isfolder(fullfile(candidates{i}, '..', '..', 'Data', 'Finance'))
        demo_dir = candidates{i};
        break
    end
end
assert(~isempty(demo_dir), "Could not locate the demo folder. " + ...
    "cd to MFTFCCA_MATLAB/demo and run again.")

addpath(genpath(fullfile(demo_dir, '..', 'src')));

data_folder = fullfile(demo_dir, '..', '..', 'Data', 'Finance');
fig_folder = string(fullfile(demo_dir, 'Figures', 'Finance'));
if ~isfolder(fig_folder), mkdir(fig_folder); end

%% Read in CPI (monthly) and take the YOY growth
% Keep only the rows tagged 'M<month>', reshaped to (year, month, value).

[value_CPI, txt_CPI] = xlsread(fullfile(data_folder, "CPI.xlsx"), 'Data');
txt_CPI = txt_CPI(2:end, :);

data_CPI = [];
for i = 1:length(value_CPI)
    mon = char(txt_CPI(i, 2));
    if mon(1) == 'M'
        data = [value_CPI(i, 1), str2double(mon(2:end)), value_CPI(i, 3)];
        data_CPI = [data_CPI; data];
    end
end
clear value_CPI txt_CPI mon data

yoy_CPI = YOY(data_CPI, 12, 1986, 2016);
yoy_CPI = reshape(yoy_CPI', [], 1);

%% Read in OIL (monthly) and take the YOY growth
% Dates arrive as 'yyyy/mm' strings.

[value_OIL, txt_OIL] = xlsread(fullfile(data_folder, "OIL.xls"), 'Data 1');
txt_OIL = txt_OIL(4:end, :);

data_OIL = [];
for i = 1:length(value_OIL)
    date = strsplit(char(txt_OIL(i,1)), '/');
    data = [str2double(string(date(1))), str2double(string(date(2))), value_OIL(i)];
    data_OIL = [data_OIL; data];
end
clear value_OIL txt_OIL

yoy_OIL = YOY(data_OIL, 12, 1986, 2016);
yoy_OIL = reshape(yoy_OIL', [], 1);

%% Read in GDP (quarterly, already a growth rate)
% Headers are 'yyyyQq', so the year and quarter are parsed off the label.

[value_GDP, txt_GPD] = xlsread(fullfile(data_folder, "GDP.xlsx"), 'Sheet1');
data_GDP = [zeros(2, length(value_GDP)); value_GDP];
for i = 1:length(value_GDP)
    date = char(txt_GPD(i));
    data_GDP(1,i) = str2double(date(1:4));
    data_GDP(2,i) = str2double(date(end));
end
clear value_GDP txt_GPD

yoy_GDP = cut_data_by_year(data_GDP', 1987, 2016);
yoy_GDP = yoy_GDP(:, 3);

%% Overview of the three YOY series

figure;
plot(yoy_CPI', '.-');
hold on
plot(1:3:30*12, yoy_GDP, '.-')
yyaxis right
plot(yoy_OIL, '.-')
hold off

%% GDP growth alone

figure("Position", [0 0 600 400]);
plot(1:3:30*12, yoy_GDP, '.-', "LineWidth", 6, "Markersize", 1, "Color", [1 0.1 0.3]);
set(gca, "fontSize", 25)
set(gca, "LineWidth", 4)
axis padded
xticks([1, 30*12]); xticklabels(["1987 Jan", "2016 Jan"])
title("YOY Growth (%)")
ylabel("GDP")
% saveas(gcf, fullfile(fig_folder, "GDP"), "epsc")

%% Time-frequency Granger causality between CPI and OIL

fs = 12;                    % monthly series
fs2 = 4;                    % quarterly series
downsample_ratio = 3;
ntrials = 1;
nobs = 30*12;

G_cause(permute(cat(3, yoy_CPI, yoy_OIL), [3,1,2]), fs, ...
    "mode", "tf", "plot", true, "report", true);
