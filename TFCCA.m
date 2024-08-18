%% TF-CCA
% 
% Estimate the mixed-frequency time-frequency CCA between X (high sampling rate) and Y (low sampling rate)
% 
% 
%% Arguments
% 
% _input_
%
%   X                       The high sampling rate data, shape (n_trials, n_high_samples)
%   Y                       The low sampling rate data, shape (n_trials, n_low_samples)
%                           X and Y should have the same number of trials and be aligned at the start point.
%   fs                      The (high) sampling rate of X
%   downsample_ratio        The ratio of the high sampling rate to the low sampling rate
%   win                     The window length to compute STFT. Number of X data points within a window.
%   time_start              The start time of the data to be analyzed, in seconds
%   time_end                The end time of the data to be analyzed, in seconds
%   real_STFT               False: Compute complex-valued STFT. True: Compute real-valued STFT.
% 
%   optional arguments         
%       calibrations        The calibration r values, shape (n_trials, 2). Default: [] 
%       lag                 The lag between X and Y, in seconds. Default: 0
%       FFT_num             The number of FFT points. Default: win
%       surrogate_num       The number of surrogates per trial. Default: int(2000/n_trials)
%       fpass               The frequency passband. Default: []
%       partial_high_set    The partial data set in high sampling rate, shape (n_trials, n_high_samples). Default: []
%       partial_low_set     The partial data set in low sampling rate, shape (n_trials, n_low_samples). Default: []
%       plot                Whether to plot the CCA Coefficients. Default: true
%       plot_interval       The interval of trials to plot. Default: int(n_trials/100)
%
% _output_
%
%   stft_f                  The frequency axis of the STFT, shape (n_freqs,)
%   r_ori                   The original TF-CCA correlations, shape (n_trials,)
%   As0                     The original TF-CCA coefficients, shape (n_freqs, n_trials)
%   r_shuffle               The surrogate TF-CCA correlations, shape (n_trials, surrogate_num)
%   As                      The surrogate TF-CCA coefficients, shape (n_freqs, n_trials * surrogate_num)
% 
%%
function [stft_f, r_ori, As0, r_shuffle, As] = TFCCA(X, Y, fs, downsample_ratio, win, time_start, time_end, real_STFT, options)

    arguments
        X (:,:) 
        Y (:,:) 
        fs (1,1)
        downsample_ratio (1,1) {mustBeInteger}
        win (1,1)
        time_start (1,1) double
        time_end (1,1) double
        real_STFT (1,1)
        options.calibrations (:,:) = []
        options.lag (1,1) double = 0
        options.FFT_num (1,1) double = win
        options.surrogate_num (1,1) {mustBeInteger} = ceil(2000/size(X,1))
        options.fpass double = []
        options.partial_high_set (:,:) = []
        options.partial_low_set (:,:) = []
        options.plot (1,1) logical = true
        options.plot_interval (1,1) {mustBeInteger} = round(size(X, 1) / 100)

    end

    ntrials = size(X, 1); x_len = size(X, 2);
    lag = options.lag;
    FFT_num = options.FFT_num;
    surrogate_num = options.surrogate_num;
    fpass = options.fpass;
    calibrations = options.calibrations;
    plot_interval = options.plot_interval;

    if isempty(calibrations)
        r_cal_max = ones(ntrials, 1);
        r_cal_min = zeros(ntrials, 1);
    else
        r_cal_max = calibrations(:, 1);
        r_cal_min = calibrations(:, 2);
    end

    partial = ~isempty(options.partial_high_set); partial_high_set = options.partial_high_set;
    partial_low = ~isempty(options.partial_low_set);
    if partial_low
        partial_low_set = options.partial_low_set;
        fs2 = fs / downsample_ratio;
        win2 = ceil(win / downsample_ratio);
        FFT_num2 = FFT_num / downsample_ratio;
    end
    assert(partial_low == 0 || partial == 0, "Only one of partial_high_set and partial_low_set can be set")

    LH_index = @(x) (x-1) * downsample_ratio + 1;
    % HL_index = @(x) ceil((x-1) / downsample_ratio) + 1;
    
    x22_r1 = ceil(time_start * fs / downsample_ratio) + 1;
    x22_r2 = floor(time_end * fs / downsample_ratio) + 1;
    lag = round(lag * fs - (win - 1) / 2);
    s1_r1 = LH_index(x22_r1) + lag;
    s1_r2 = LH_index(x22_r2) + lag;

    assert(s1_r1 >= 1, "time_start is too small for the lag")
    assert(x_len - win + 1 >= s1_r2, "time_end is too large for the lag")


    % -------- Original --------

%     clf;
    As0 = [];
    r_ori = zeros(1, ntrials);

    for i = progress(1:ntrials)
        
        x_high = X(i,1:end); x_low = Y(i,1:end);

        [s1, stft_f, ~] = STFT(x_high, fs, win, win-1, FFT_num, 'plot', 0, "fpass", fpass);
        s1 = zscore(s1')';
        if real_STFT
            s1 = abs(s1);
        end        
        s1 = s1(:, s1_r1 : downsample_ratio : s1_r2);

        if partial
            x_partial = partial_high_set(i,:);       
            s3 = STFT(x_partial, fs, win, win-1, FFT_num, 'plot', 0, "fpass", fpass);
            s3 = zscore(s3')';
            if real_STFT
                s3 = abs(s3);
            end
            s3 = s3(:, s1_r1 : downsample_ratio : s1_r2);
            [A,~,ri] = pCCA((s1'), x_low(x22_r1:x22_r2)', (s3'));
        elseif partial_low
            x_partial = partial_low_set(i,:);
            s3 = STFT(x_partial, fs2, win2, win2-1, FFT_num2, 'plot', 0, "fpass", fpass);
            s3 = zscore(s3')';
            if real_STFT
                s3 = abs(s3);
            end
            s3 = s3(:, HL_index(s1_r1) : downsample_ratio : HL_index(s1_r2));
            [A,~,ri] = pCCA((s1'), x_low(x22_r1:x22_r2)', (s3'));
        else
            [A,~,ri] = canoncorr(s1',(x_low(x22_r1:x22_r2)'));
        end
        
        ri = (ri - r_cal_min(i)) / (r_cal_max(i) - r_cal_min(i));
        A = A/mean(abs(A));
        
        As0 = [As0, A];
        r_ori(i) = ri;

    end

    if options.plot
        plot_CCA_coeff(stft_f, As0, "plot_interval", plot_interval)
    end



    % -------- Surrogate --------


    As = [];
    r = zeros(ntrials, 2);
    r_shuffle = zeros(ntrials, surrogate_num);
    for i = progress(1:ntrials)
        x_high = X(i, 1:end); x_low = Y(i, 1:end); 
        x_high_set = phaseran(x_high', surrogate_num);
        r_j = zeros(1,surrogate_num);

        for j = 1:surrogate_num

            x_high = x_high_set(:,:,j)';
            x_low = phaseran(x_low', 1)';
            s1 = STFT(x_high, fs, win, win-1, FFT_num, 'plot',0, "fpass", fpass);
            s1 = zscore(s1')';
            if real_STFT
                s1 = abs(s1);
            end            
            s1 = s1(:, s1_r1 : downsample_ratio : s1_r2);

            if partial
                x_partial = partial_high_set(i,:);
                s3 = STFT(x_partial, fs, win, win-downsample_ratio, FFT_num, 'plot', 0, "fpass", fpass);
                s3 = zscore(s3')';
                if real_STFT
                    s3 = abs(s3);
                end
                s3 = s3(:, s1_r1 : downsample_ratio : s1_r2);
                [A,~,ri] = pCCA(s1', x_low(x22_r1:x22_r2)', (s3'));
            elseif partial_low
                x_partial = partial_low_set(i,:);
                s3 = STFT(x_partial, fs2, win2, win2-1, FFT_num2, 'plot', 0, "fpass", fpass);
                s3 = zscore(s3')';
                if real_STFT
                    s3 = abs(s3);
                end
                s3 = s3(:, HL_index(s1_r1) : downsample_ratio : HL_index(s1_r2));
                [A,~,ri] = pCCA(s1', x_low(x22_r1:x22_r2)', (s3'));
            else
                [A,~,ri] = canoncorr(s1',x_low(x22_r1:x22_r2)');
            end
            ri = (ri - r_cal_min(i)) / (r_cal_max(i) - r_cal_min(i));
            A = A / mean(abs(A));
            As = [As, A];
            
            r_j(j) = ri;
        end

        r(i,1) = mean(r_j);
        r(i,2) = var(r_j);

        r_shuffle(i,:) = r_j;

    end

    % -------- Significance test --------

    [h, p] = ttest((r_ori'-r(:,1))./sqrt(r(:,2)));
    if h == 1
        fprintf("ttest: Causal, p=%.2f", p)
    else
        fprintf("ttest: Non-causal, p=%.2f", p)
    end

    [h, p] = kstest2(r_ori, r(:,1));
    if h == 1
        fprintf("KStest: Causal, p=%.2f", p)
    else
        fprintf("KStest: Non-causal, p=%.2f", p)
    end

end