function MF_xcorr(X, Y, fs, downsample_ratio, maxlag1, maxlag2)

maxlag1 = maxlag1 * fs / downsample_ratio; maxlag2 = maxlag2 * fs / downsample_ratio;

ntrials = size(X, 1);
x22_r1 = maxlag1 + 1;
x22_r2 = floor(size(X, 2) / downsample_ratio) - maxlag2;
x1_r1 = maxlag1 * downsample_ratio + 1;
x1_r2 = (x22_r2 - 1) * downsample_ratio + 1;

maxlag1 = maxlag1 * downsample_ratio; maxlag2 = maxlag2 * downsample_ratio;
lags = -maxlag1 : 1 : maxlag2;
r = zeros(ntrials, length(lags));

for i = 1:ntrials
    
    x_high = X(i,:); x_low = Y(i,x22_r1:x22_r2);
    
    for l_index = 1:length(lags)
        lag = lags(l_index);
        
        ri = corrcoef(x_high(x1_r1+lag : downsample_ratio : x1_r2+lag)', x_low');
        r(i, l_index) = ri(1,2);           
    end
        
end

lags = lags ./ fs;
figure('Position', [100, 100, 600, 400]); hold on; set(gca, "fontSize", 18);

if ntrials > 1
    r = [mean(r);sqrt(var(r))/ sqrt(ntrials)]; 
    fill([lags lags(end:-1:1)], [r(1,:)+r(2,:), r(1,end:-1:1)-r(2,end:-1:1)], 'r', 'FaceColor', [0.4660 0.6740 0.1880] +0.3, 'EdgeColor','none')
end

plot(lags, r(1,:), 'color', '#77AC30', "LineWidth", 6)

hold off
axis padded
xlabel("Time lag (s)");
ylabel("XCorr")
set(gca, "fontSize", 25)
set(gca, "LineWidth", 4)

end