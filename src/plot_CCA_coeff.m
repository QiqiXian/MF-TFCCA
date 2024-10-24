function fig = plot_CCA_coeff(stft_f1, As1, options)
        
arguments
    stft_f1
    As1
    options.As2 = nan
    options.plot_interval (1,1) {mustBeInteger} = round(size(As1, 2) / 100)
end

ntrials = size(As1, 2);
plot_interval = options.plot_interval;

if isnan(options.As2)
    
    figure("Position", [0,0,500,400]);
    color = colormap("summer(" + num2str(ntrials) + ")"); % summer, hsv
    ax1 = gca; hold(ax1, "on"); 

    for i = 1:plot_interval:ntrials
        plot(ax1, abs(stft_f1), abs(As1(:, i)), LineWidth=1.5, Color=[color(i,:),0.3]);
    end

    plot(ax1, stft_f1, (mean(abs(As1), 2) + std(abs(As1), 0, 2)), '--', 'Color', "#228B22", 'LineWidth', 2, "MarkerSize", 15);
    plot(ax1, stft_f1, (mean(abs(As1), 2) - std(abs(As1), 0, 2)), '--', 'Color', "#228B22", 'LineWidth', 2, "MarkerSize", 15);
    fill([stft_f1; flip(stft_f1)], [mean(abs(As1), 2) + std(abs(As1), 0, 2); flip(mean(abs(As1), 2) - std(abs(As1), 0, 2))], "k", 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    h0 = plot(ax1, stft_f1, (mean(abs(As1), 2)), '-', 'Color', "#228B22", 'LineWidth', 6, "MarkerSize", 15);

    set(ax1, "fontSize", 24); set(ax1, "LineWidth", 4); 
    xlabel(ax1, 'Freq(Hz)', fontsize=30)
    ylabel(ax1, 'CCA Coeff. (abs)', fontsize=30)
    box(ax1, 'on'); hold(ax1, 'off');
        
else
        
    As2 = options.As2;
    figure("Position", [0,0,800,400]);
    color = colormap(sprintf("summer(%d)", ntrials)); % summer, hsv
    color2 = colormap(sprintf("autumn(%d)", ntrials));
    ax1 = subplot(121);
    ax2 = subplot(122);
    hold(ax1, "on"); hold(ax2, "on");


    for i = 1:plot_interval:ntrials
        plot(ax1, abs(stft_f1), abs(As1(:, i)), LineWidth=1.5, Color=[color(i,:),0.3]);
        plot(ax2, abs(stft_f1), abs(As2(:, i)), LineWidth=1.5, Color=[color2(i,:),0.3]);
    end

    h0 = plot(ax1, stft_f1, (mean(abs(As1), 2)), '-', 'Color', "#228B22", 'LineWidth', 6, "MarkerSize", 15);
    hh = plot(ax2, stft_f1, (mean(abs(As2), 2)), '-', 'Color', "#D95319", 'LineWidth', 6, "MarkerSize", 15);
    % legend([h0, hh], ["X to Y (n="+num2str(ntrials)+")", "Y to X (n="+num2str(ntrials)+")"], "Location", "best")

    set(ax1, "fontSize", 24); set(ax1, "LineWidth", 4); set(ax2, "fontSize", 24); set(ax2, "LineWidth", 4); 
    xlabel(ax1, 'Freq(Hz)', fontsize=30)
    ylabel(ax1, 'CCA Coeff. (abs)', fontsize=30)
    title(ax1, "X to Y"); title(ax2, "Y to X")
    box(ax1, 'on'); box(ax2, 'on'); hold(ax1, 'off');
end

fig = gcf;

end