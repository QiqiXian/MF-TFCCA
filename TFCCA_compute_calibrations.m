%% Compute correlation calibrations from X

function [calibrations] = TFCCA_compute_calibrations(X, fs, downsample_ratio, win, time_start, time_end, maxlag1, maxlag2, real_STFT, options)

    arguments
        X (:,:) 
        fs (1,1)
        downsample_ratio (1,1)
        win (1,1)
        time_start (1,1) double
        time_end (1,1) double
        maxlag1 (1,1) double
        maxlag2 (1,1) double
        real_STFT (1,1) logical
        options.FFT_num (1,1) double = win
        options.fpass = []
        options.plot (1,1) logical = true
    end

    fpass = options.fpass;
    FFT_num = options.FFT_num;
    
    lag_max = ceil(maxlag1 * fs) + ceil((win - 1) / 2 );
    lag_max2 = floor(maxlag2 * fs) - ceil((win - 1) / 2 );
    
    x22_r1 = ceil(time_start * fs / downsample_ratio) + 1;
    x22_r2 = floor(time_end * fs / downsample_ratio) + 1;
    LH_index = @(x) (x-1) * downsample_ratio + 1;
    s1_r1 = LH_index(x22_r1);
    s1_r2 = LH_index(x22_r2);

    x_len = size(X,2);
    assert (s1_r1 > 0, "time_start is too small")
    assert (s1_r2 <= x_len, "time_end is too large")
    if s1_r1 < lag_max
        lag_max = s1_r1 - 1;
        fprintf("maxlag1 is too large, set to %d", lag_max)
    end
    % assert (s1_r1 >= lag_max, "maxlag1 is too large")
    sx_len = x_len - win + 1;
    if sx_len - s1_r2 < lag_max2
        lag_max2 = sx_len - s1_r2;
        fprintf("maxlag2 is too large, set to %d", lag_max2)
    end
    % assert(sx_len - s1_r2 >= lag_max2, "maxlag2 is too large")
    
    lags = (-lag_max + 1 : lag_max2) / fs  + (win - 1) / 2 / fs;
    r = zeros(size(X,1), lag_max+lag_max2);
    
    X_downsampled = X(:, 1:downsample_ratio:end);

    for i = 1:size(X,1)
        x_high = X(i,:); x_low = X_downsampled(i,:);
    
        [s1, ~, ~] = STFT(x_high, fs, win, win-1, FFT_num, 'plot',0, "fpass", fpass);    
        s1 = zscore(s1')';
        if real_STFT
            s1 = abs(s1);
        end
    
        for lag = -lag_max + 1 : lag_max2
            [~,~,ri] = canoncorr((s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)'),x_low(x22_r1:x22_r2)');
            r(i, lag_max + lag) = ri;
        end
        
    end
        
    r_cal_max = max(r, [], 2);
    r_cal_min = min(r, [], 2);
    calibrations = [r_cal_max, r_cal_min];

    r = (r - r_cal_min) ./ (r_cal_max - r_cal_min);
    
    if options.plot
        figure('Position', [100, 100, 600, 400]); hold on; set(gca, "fontSize", 18);
        if size(X,1) > 1
            r = [mean(r);sqrt(var(r)) ];
            fill([lags lags(end:-1:1)], [r(1,:)+r(2,:), r(1,end:-1:1)-r(2,end:-1:1)], 'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor','none')
        end
        yLimits = ylim(); plot([-win/fs/2, -win/fs/2], yLimits, "k--", "LineWidth", 3); plot([win/fs/2, win/fs/2], yLimits, "k--", "LineWidth", 3);
        plot(lags, r(1,:), "LineWidth", 8)
        hold off; axis padded; set(gca, "fontSize", 25); set(gca, "LineWidth", 4)
        xlabel("Time lag (s)", fontsize=30); ylabel("Can. Corr.", fontsize=30)
    end
        

end