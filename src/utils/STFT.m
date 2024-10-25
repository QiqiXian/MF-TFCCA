function [s,f,t] = STFT(x, fs, window_length, OverlapLength, FFTL, options)

arguments
    x
    fs
    window_length
    OverlapLength
    FFTL
    options.plot (1,1) = true
    options.fpass = []
end

% windowl < FFTL < xlength

[s,f,t] = stft(x, fs, Window=hann(window_length), OverlapLength=OverlapLength,FFTLength=FFTL);
if isempty(options.fpass)
    fpass = [0 fs/2];
else
    fpass = options.fpass;
end
for i = 1:length(f)
    if f(i) >= fpass(1) - 1e-6
        f_start = i;
        break
    end
end
for j = i : length(f)
    if (f(j) >= fpass(2) - 1e-6) || (j==length(f))
        f_end = j;
        break
    end
end
f = f(f_start:f_end); s = s(f_start:f_end, :);
if options.plot
    cla; clf; figure;
    ax = gca;
    sdb = mag2db(abs(s));
    imagesc(t, f(end:-1:1), sdb(end:-1:1,:))
%     cc = max(sdb(:))+[-60 0];
%     ax.CLim = cc;
%     ax.YDir = 'reverse';
    axis xy
    view(2)
    colormap('turbo')
    colorbar
    xlabel('time(s)')
    ylabel('freq(Hz)')
    title('STFT spectrum (db)')
    set(ax, "FontSize", 20)
end

end
