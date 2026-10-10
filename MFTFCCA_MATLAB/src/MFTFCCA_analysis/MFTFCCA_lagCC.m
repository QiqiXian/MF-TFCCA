%% TFCCA_lagCC
% Compute the TFCCA lag-CC profile
% 
%%

function [lags, r_0, r_su] = MFTFCCA_lagCC(X, Y, fs, downsample_ratio, win, time_start, time_end, maxlag1, maxlag2, real_STFT, options)

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
        options.calibrations (:,:) = []
        options.FFT_num (1,1) double = win
        options.surrogate_num (1,1) {mustBeInteger} = round(2000/size(X,1))
        options.fpass double = []
        options.partial_high_set = []
        options.partial_low_set = []
        options.partial_lag_set = []
        options.n_lags (1,1) double = NaN
        options.plot (1,1) logical = true
        options.fig = []
    end

    ntrials = size(X, 1); x_len = size(X, 2);
    FFT_num = options.FFT_num;
    surrogate_num = options.surrogate_num;
    fpass = options.fpass;
    calibrations = options.calibrations;

    if isscalar(calibrations) && calibrations == 1
        disp("Compute calibrations...")
        calibrations = MFTFCCA_compute_calibrations(X, fs, downsample_ratio, ...
            win, time_start, time_end, maxlag1, maxlag2, real_STFT, ...
            "FFT_num", FFT_num, "fpass", fpass, "plot", 0);
        disp("Calibrations computed")
    end
    if isempty(calibrations) || (isscalar(calibrations) && calibrations == 0)
        r_cal_max = ones(ntrials, 1);
        r_cal_min = zeros(ntrials, 1);
    else
        r_cal_max = calibrations(:, 1);
        r_cal_min = calibrations(:, 2);
    end

    partial = ~isempty(options.partial_high_set); 
    if partial
        partial_high_set = convert_to_3d_array(options.partial_high_set);
    end

    fs2 = fs / downsample_ratio;
    win2 = ceil(win / downsample_ratio);
    FFT_num2 = FFT_num / downsample_ratio;

    partial_low = ~isempty(options.partial_low_set);
    if partial_low
        partial_low_set = convert_to_3d_array(options.partial_low_set);
    end

    % Raw low-rate signals whose recent past is regressed out at every
    % analysis point: at point t of the Y window, the taps
    % z(t-1) ... z(t-n_lags) are partialled out. The taps are indexed by Y's
    % own time base, so the block does not move with the lag and z(t) itself
    % is never included. Equivalent to removing any FIR-filtered
    % contribution of z over that span, without going through its spectrum.
    partial_lag = ~isempty(options.partial_lag_set);
    n_lags = options.n_lags;
    if isnan(n_lags)
        n_lags = win2;
    end
    if n_lags == 0
        partial_lag = 0;
    end
    if partial_lag
        partial_lag_set = convert_to_3d_array(options.partial_lag_set);
        assert(n_lags >= 1, "n_lags must be at least 1")
        assert(size(partial_lag_set, 3) >= floor(x_len / downsample_ratio), ...
            "partial_lag_set must be on the low-rate grid (N_trials x time/downsample_ratio)")
    end

    % assert(partial_low == 0 || partial == 0, "Only one of partial_high_set and partial_low_set can be set")

    LH_index = @(x) (x-1) * downsample_ratio + 1;
    HL_index = @(x) ceil((x-1) / downsample_ratio) + 1;
   
    x22_r1 = ceil(time_start * fs / downsample_ratio) + 1;
    x22_r2 = floor(time_end * fs / downsample_ratio) + 1;
    s1_r1 = LH_index(x22_r1);
    s1_r2 = LH_index(x22_r2);

    assert (x22_r2 <= x_len / downsample_ratio, "time_end is too large")
    assert (x22_r1 <= x22_r2, "time_start is larger than time_end")

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
    
    if partial_lag
        assert(x22_r1 - n_lags >= 1, ...
            "Not enough samples before time_start for %d lag taps: increase " + ...
            "time_start or reduce n_lags.", n_lags)
        assert(x22_r2 <= size(partial_lag_set, 3), ...
            "partial_lag_set is shorter than the Y analysis window.")
    end

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
        
        % % Generate partial dataset

        
        if partial
            s3_high = [];
            for k = 1:size(partial_high_set, 1)
                x_partial = squeeze(partial_high_set(k,:,:)); 
                x_partial = x_partial(i,:);
                [s3_, ~, ~] = STFT(x_partial, fs, win, win-1, FFT_num, "plot", 0, "fpass", fpass);
                s3_ = zscore(s3_')';
                if real_STFT
                    s3_ = abs(s3_);
                end 
                s3_high = [s3_high; s3_];
            end
        end
        if partial_low
            s3_low = [];
            for k = 1:size(partial_low_set, 1)
                x_partial = squeeze(partial_low_set(k,:,:)); 
                x_partial = x_partial(i,:);
                [s3_, ~, ~] = STFT(x_partial, fs2, win2, win2-1, FFT_num2, "plot", 0, "fpass", fpass);
                s3_ = zscore(s3_')';
                if real_STFT
                    s3_ = abs(s3_);
                end 
                s3_low = [s3_low; s3_];
            end
            % s3 = s3(abs(mean(abs(s3),2)) > 5e-2, :); % remove low power freqs
        end


        if partial_lag
            z_lag = reshape(partial_lag_set(:, i, :), size(partial_lag_set, 1), []);
            s3_lag = build_lag_block(z_lag, x22_r1, x22_r2, n_lags);
        end

        for lag = -lag_max + 1 : lag_max2
            s3 = [];
            if partial
                s3 = [s3; s3_high(:, s1_r1+lag: downsample_ratio : s1_r2+lag)];
            end
            if partial_low
                s3 = [s3; s3_low(:, HL_index(s1_r1+lag):HL_index(s1_r2+lag))];
            end
            if partial_lag
                s3 = [s3; s3_lag];
            end

            if partial || partial_low || partial_lag
                [~,~,ri] = pCCA(s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)', x_low(x22_r1:x22_r2)', s3');
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
            s3_high = [];
            for k = 1:size(partial_high_set, 1)
                x_partial = squeeze(partial_high_set(k,:,:)); 
                x_partial = x_partial(i,:);
                [s3_, ~, ~] = STFT(x_partial, fs, win, win-1, FFT_num, "plot", 0, "fpass", fpass);
                s3_ = zscore(s3_')';
                if real_STFT
                    s3_ = abs(s3_);
                end 
                s3_high = [s3_high; s3_];
            end
        end

        if partial_low
            s3_low = [];
            for k = 1:size(partial_low_set, 1)
                x_partial = squeeze(partial_low_set(k,:,:)); 
                x_partial = x_partial(i,:);
                [s3_, ~, ~] = STFT(x_partial, fs2, win2, win2-1, FFT_num2, "plot", 0, "fpass", fpass);
                s3_ = zscore(s3_')';
                if real_STFT
                    s3_ = abs(s3_);
                end
                s3_low = [s3_low; s3_];
            end
            % s3 = s3(abs(mean(abs(s3),2)) > 5e-2, :); % remove low power freqs
        end

        if partial_lag
            z_lag = reshape(partial_lag_set(:, i, :), size(partial_lag_set, 1), []);
            s3_lag = build_lag_block(z_lag, x22_r1, x22_r2, n_lags);
        end

        for s = 1:surrogate_num
            x_high = x_high_su(:,:,s)';
            % x_low = phaseran(x_low', 1)';

            [s1, ~, ~] = STFT(x_high, fs, win, win-1, FFT_num, 'plot',0, "fpass", fpass);
        
            s1 = zscore(s1')';
            if real_STFT
                s1 = abs(s1);
            end
            
            for lag = -lag_max + 1 : lag_max2
                s3 = [];
                if partial
                    s3 = [s3; s3_high(:, s1_r1+lag: downsample_ratio : s1_r2+lag)];
                end
                if partial_low
                    s3 = [s3; s3_low(:, HL_index(s1_r1+lag):HL_index(s1_r2+lag))];
                end
                if partial_lag
                    s3 = [s3; s3_lag];
                end

                if partial || partial_low || partial_lag
                    [~,~,ri] = pCCA(s1(:, s1_r1+lag: downsample_ratio : s1_r2+lag)', x_low(x22_r1:x22_r2)', s3');
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
        % determine target axes: if options.fig is an axes handle, use it; else create figure
        ax = [];
        if isfield(options, "fig") && isa(options.fig, 'matlab.graphics.axis.Axes') && isvalid(options.fig)
            ax = options.fig;
            hold(ax, "on");
        else
            fig = figure('Position', [100, 100, 800, 400]);
            ax = axes(fig); hold(ax, "on");
        end

        norm_factor = sqrt(ntrials);
        % norm_factor = 1;

        r = r_su;
        fill(ax, [lags lags(end:-1:1)], [r(1,:)+r(2,:)/norm_factor, r(1,end:-1:1)-r(2,end:-1:1)/norm_factor], ...
            'b', 'FaceColor', [0.8 0.8 1], 'EdgeColor','none');

        if ntrials > 1
            r = r_0;
            fill(ax, [lags lags(end:-1:1)], [r(1,:)+r(2,:)/norm_factor , r(1,end:-1:1)-r(2,end:-1:1)/norm_factor], ...
                'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor','none');
        end

        colororder(ax, [0    0.4470    0.7410; ...
            0.8500    0.3250    0.0980;])
        for r_ = {r_su, r_0}
            r = cell2mat(r_);
            plot(ax, lags, r(1,:), "LineWidth", 6);
        end

        yLimits = ylim(ax);

        plot(ax, [-win/fs/2, -win/fs/2], yLimits, "k--", "LineWidth", 3);
        plot(ax, [win/fs/2, win/fs/2], yLimits, "k--", "LineWidth", 3);

        hold(ax, "off")
        % legend(ax, ["", "", "Surrogate", "Original"], "Location", "bestoutside");
        axis(ax, "padded")
        xlabel(ax, "Time lag (s)");
        ylabel(ax, "Can. Corr")
        set(ax, "fontSize", 25, "LineWidth", 4)

    end

end


function blk = build_lag_block(Zi, t1, t2, n_lags)
% Past-sample design block for the raw partial-out signals.
%
%   Zi     - K x time matrix, one row per partial-out variable (one trial)
%   t1, t2 - first and last point of the Y analysis window (low-rate indices)
%   n_lags - number of taps
%
% Row (k-1)*n_lags + m holds variable k delayed by m samples, i.e. the
% column for analysis point t carries Zi(k, t-m). Taps run from t-1 back to
% t-n_lags; the current sample z(t) is not included.

    K = size(Zi, 1);
    blk = zeros(K * n_lags, t2 - t1 + 1);
    row = 0;
    for k = 1:K
        for m = 1:n_lags
            row = row + 1;
            blk(row, :) = Zi(k, t1-m : t2-m);
        end
    end
end
