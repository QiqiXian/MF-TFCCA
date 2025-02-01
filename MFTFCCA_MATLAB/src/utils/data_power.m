function [S, f] = data_power(data, params)

figure("Position", [0,0,300,250]); hold on;

if (params.trialave == 1)
    
    [S,f, serr]=mtspectrumc(squeeze(data), params);
%     figure;

    % plot(f,serr, 'k-')
    fill([f f(end:-1:1)], mag2db([serr(1,:), serr(2,end:-1:1)]), 'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor','none')
    
else
    [S,f]=mtspectrumc(squeeze(data), params);
    
end

% SS = smoothdata(mag2db(S), "gaussian", 160);
SS = mag2db(S);
h = plot(f(1:end), SS(1:end), 'k-', "LineWidth", 4);
% set (h, 'linesmoothing', 'on');
% fill([f f(end:-1:1)], ([serr(1,:), serr(2,end:-1:1)]), 'r', 'FaceColor', [1 0.8 0.8], 'EdgeColor','none')
% plot(f, (S), 'k-', "LineWidth", 0.3)
% hold off

xlabel('Freq (Hz)')
ylabel('PSD (dB)')
title('Power spectrum')
end