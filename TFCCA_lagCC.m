%% TFCCA_lagCC
% Compute the TFCCA lag-CC profile
% 
%%

function [lags, r_0, r_su] = TFCCA_lagCC(X, Y, fs, downsample_ratio, win, time_start, time_end, maxlag1, maxlag2, real_STFT, calibrations, options)

    arguments
        X (:,:) 
        Y (:,:) 
        fs (1,1)
        downsample_ratio (1,1)
        win (1,1)
        time_start (1,1) double
        time_end (1,1) double
        maxlag1 (1,1) double
        maxlag2 (1,1) double
        real_STFT (1,1) logical
        calibrations (:,:) = []
        options.FFT_num (1,1) double = win
        options.surrogate_num (1,1) {mustBeInteger} = round(2000/size(X,1))
        options.fpass double = []
        options.partial_high_set (:,:) = []
        options.partial_low_set (:,:) = []
        options.plot (1,1) logical = true
    end

    ntrials = size(X, 1); x_len = size(X, 2);
    FFT_num = options.FFT_num;
    surrogate_num = options.surrogate_num;
    fpass = options.fpass;

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
    HL_index = @(x) ceil((x-1) / downsample_ratio) + 1;
    
    x22_r1 = ceil(time_start * fs / downsample_ratio) + 1;
    x22_r2 = floor(time_end * fs / downsample_ratio) + 1;
    s1_r1 = LH_index(x22_r1);
    s1_r2 = LH_index(x22_r2);

    lag_max = ceil(maxlag1 * fs) + ceil((win - 1) / 2 );
    lag_max2 = floor(maxlag2 * fs) - ceil((win - 1) / 2 );
    
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
    % assert(size(s1, 2) - s1_r2 >= lag_max2)
    
    lags = (-lag_max + 1 : lag_max2) / fs  + (win - 1) / 2 / fs;
    

    % -------- Original --------

    r = zeros(ntrials, lag_max+lag_max2);

    for i = 1:ntrials
      
        x_high = X(i,:); x_low = Y(i,:);

        [s1, ~, ~] = STFT(x_high, fs, win, win-1, FFT_num, 'plot',0, "fpass", fpass);    
        s1 = zscore(s1')';
        if real_STFT
            s1 = abs(s1);
        end
        
        if partial
            x_partial = partial_high_set(i,:);
            [s3, ~, ~] = STFT(x_partial, fs, win, win-1, FFT_num, "plot", 0, "fpass", fpass);
            s3 = zscore(s3')';
            if real_STFT
                s3 = abs(s3);
            end 

        elseif partial_low
            x_partial = partial_low_set(i,:);
            [s3, ~, ~] = STFT(x_partial, fs2, win2, win2-1, FFT_num2, "plot", 0, "fpass", fpass);
            s3 = zscore(s3')';
            if real_STFT
                s3 = abs(s3);
            end 
            % s3 = s3(abs(mean(abs(s3),2)) > 5e-2, :); % remove low power freqs
        end

        for lag = -lag_max + 1 : lag_max2

            if partial
                [~,~,ri] = pCCA(s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)', x_low(x22_r1:x22_r2)', s3(:, s1_r1+lag: downsample_ratio : s1_r2+lag)');
            elseif partial_low
                [~,~,ri] = pCCA(s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)', x_low(x22_r1:x22_r2)', s3(:, HL_index(s1_r1+lag):HL_index(s1_r2+lag))');
            else
                [~,~,ri] = canoncorr((s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)'),x_low(x22_r1:x22_r2)');
            end
            r(i, lag_max + lag) = ri;
        end
            
    end

    r = (r - r_cal_min(1:ntrials)) ./ (r_cal_max(1:ntrials) - r_cal_min(1:ntrials));
    if ntrials > 1
        r = [mean(r);sqrt(var(r)) ]; %  / sqrt(ntrials)
    end

    r_0 = r;

    
    % -------- Surrogate --------

    r = zeros(ntrials, lag_max+lag_max2, surrogate_num);

    for i = 1:ntrials
        
        x_high = X(i,:); x_low = Y(i,:);
        x_high_su = phaseran(x_high', surrogate_num);
        
        if partial
            x_partial = partial_high_set(i,:);
            [s3, ~, ~] = STFT(x_partial, fs, win, win-1, FFT_num, "plot", 0, "fpass", fpass);
            s3 = zscore(s3')';
            if real_STFT
                s3 = abs(s3);
            end 

        elseif partial_low
            x_partial = partial_low_set(i,:);
            [s3, ~, ~] = STFT(x_partial, fs2, win2, win2-1, FFT_num2, "plot", 0, "fpass", fpass);
            s3 = zscore(s3')';
            if real_STFT
                s3 = abs(s3);
            end 
            % s3 = s3(abs(mean(abs(s3),2)) > 5e-2, :); % remove low power freqs
        end
        
        for s = 1:surrogate_num
            x_high = x_high_su(:,:,s)';
            x_low = phaseran(x_low', 1)';

            [s1, ~, ~] = STFT(x_high, fs, win, win-1, FFT_num, 'plot',0, "fpass", fpass);
        
            s1 = zscore(s1')';
            if real_STFT
                s1 = abs(s1);
            end
            
            for lag = -lag_max + 1 : lag_max2
                if partial
                    [~,~,ri] = pCCA(s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)', x_low(x22_r1:x22_r2)', s3(:, s1_r1+lag: downsample_ratio : s1_r2+lag)');
                elseif partial_low
                    [~,~,ri] = pCCA(s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)', x_low(x22_r1:x22_r2)', s3(HL_index(s1_r1+lag):HL_index(s1_r2+lag))');
                else
                    [~,~,ri] = canoncorr((s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)'),x_low(x22_r1:x22_r2)');
                end
                
                r(i, lag_max + lag, s) = ri;
            end
        end
    end

    r = (r - r_cal_min) ./ (r_cal_max - r_cal_min);
    % r_su_all = r;
    r = [mean(r, [1,3]);sqrt(var(r, [], [1,3])) ]; %  / sqrt(ntrials)
    r_su = r;
    

    % -------- Plot --------

    if options.plot

        figure('Position', [100, 100, 800, 400]); hold on;

        norm_factor = sqrt(ntrials);

        r = r_su;
        fill([lags lags(end:-1:1)], [r(1,:)+r(2,:)/norm_factor, r(1,end:-1:1)-r(2,end:-1:1)/norm_factor], ...
        'b', 'FaceColor', [0.8 0.8 1], ...
            'EdgeColor','none')

        if ntrials > 1
            r = r_0;
            fill([lags lags(end:-1:1)], [r(1,:)+r(2,:)/norm_factor , r(1,end:-1:1)-r(2,end:-1:1)/norm_factor], ...
            'r', 'FaceColor', [1 0.8 0.8], ...
                'EdgeColor','none')
        end
        
        colororder([0    0.4470    0.7410; ...
            0.8500    0.3250    0.0980;])
        for r_ = {r_su, r_0}
            r = cell2mat(r_);
            plot(lags, r(1,:), "LineWidth", 6);
        end

        yLimits = ylim();

        plot([-win/fs/2, -win/fs/2], yLimits, "k--", "LineWidth", 3);
        plot([win/fs/2, win/fs/2], yLimits, "k--", "LineWidth", 3);

        hold off
        legend(["", "", "Surrogate", "Original"], "Location", "bestoutside");
        axis padded
        xlabel("Time lag (s)");
        ylabel("Can. Corr")
        set(gca, "fontSize", 25)
        set(gca, "LineWidth", 4)

    end

end