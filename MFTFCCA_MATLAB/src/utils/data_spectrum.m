
function [S, t, f] = data_spectrum(data, params, win, options)

arguments
    data
    params
    win
    options.plot (1,1) = true
end

[S, t, f] = mtspecgramc(squeeze(data), win, params);
if options.plot
    figure;
    ax = gca;
    sdb = mag2db(abs(S));
    imagesc(t, f, sdb'); colormap('turbo');
%     cc = max(sdb(:))+[-100 0]; ax.CLim = cc; 
    view(2); colorbar;
    axis xy
    xlabel('time(s)')
    ylabel('freq(Hz)')
    title('spectrum (db)')
end

end