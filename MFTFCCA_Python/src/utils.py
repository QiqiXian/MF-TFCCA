import numpy as np
import matplotlib.pyplot as plt

def MF_xcorr(X, Y, fs, downsample_ratio, maxlag1, maxlag2):
    """
    Parameters:
        X (ndarray): High-frequency signal (ntrials x time).
        Y (ndarray): Low-frequency signal (ntrials x time).
        fs (float): Sampling frequency.
        downsample_ratio (int): Downsampling ratio.
        maxlag1 (float): Maximum lag 1 (in seconds).
        maxlag2 (float): Maximum lag 2 (in seconds).
    """
    # Convert lags to samples
    maxlag1 = int(maxlag1 * fs / downsample_ratio)
    maxlag2 = int(maxlag2 * fs / downsample_ratio)
    
    ntrials = X.shape[0]
    x22_r1 = maxlag1 + 1
    x22_r2 = int(np.floor(X.shape[1] / downsample_ratio)) - maxlag2
    x1_r1 = maxlag1 * downsample_ratio + 1
    x1_r2 = (x22_r2 - 1) * downsample_ratio + 1
    
    maxlag1 = maxlag1 * downsample_ratio
    maxlag2 = maxlag2 * downsample_ratio
    lags = np.arange(-maxlag1, maxlag2 + 1)
    r = np.zeros((ntrials, len(lags)))
    
    # Compute cross-correlation for each trial
    for i in range(ntrials):
        x_high = X[i, :]
        x_low = Y[i, x22_r1 - 1:x22_r2]  # Python uses 0-based indexing
        
        for l_index, lag in enumerate(lags):
            x_high_lagged = x_high[x1_r1 - 1 + lag::downsample_ratio][:x22_r2 - x22_r1 + 1]
            ri = np.corrcoef(x_high_lagged, x_low)[0, 1]
            r[i, l_index] = ri
    
    # Convert lags to seconds
    lags = lags / fs
    
    # Plot results
    plt.figure(figsize=(8, 6))
    plt.rcParams.update({'font.size': 18})
    
    if ntrials > 1:
        r_mean = np.mean(r, axis=0)
        r_std = np.std(r, axis=0) / np.sqrt(ntrials)
        plt.fill_between(lags, r_mean - r_std, r_mean + r_std, color=[0.4660, 0.6740, 0.1880, 0.3], edgecolor='none')
    else:
        r_mean = r[0, :]
    
    plt.plot(lags, r_mean, color='#77AC30', linewidth=6)
    plt.xlabel("Time lag (s)")
    plt.ylabel("XCorr")
    plt.gca().tick_params(axis='both', labelsize=25)
    plt.gca().spines['bottom'].set_linewidth(4)
    plt.gca().spines['left'].set_linewidth(4)
    plt.grid(True)
    plt.show()