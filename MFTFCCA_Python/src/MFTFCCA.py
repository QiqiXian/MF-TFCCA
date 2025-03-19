import numpy as np
import matplotlib.pyplot as plt
from sklearn.cross_decomposition import CCA
from scipy.signal import stft, hann
from matplotlib.colors import Normalize
import numpy as np
from scipy.stats import ttest_ind, ks_2samp

def canoncorr(X, Y):
    """
    Parameters:
        X (ndarray): High-frequency signal (n_samples x n_features).
        Y (ndarray): Low-frequency signal (n_samples x n_features).
    
    Returns:
        U (ndarray): Canonical scores for X.
        V (ndarray): Canonical scores for Y.
        r (float): Canonical correlation coefficient.
    """
    cca = CCA(n_components=1)
    cca.fit(X, Y)
    U, V = cca.transform(X, Y)
    r = np.corrcoef(U.T, V.T)[0, 1]
    return U, V, r

def zscore(x):
    return (x - np.mean(x, axis=0)) / np.std(x, axis=0)

def STFT(x, fs, window_length, OverlapLength, FFTL, plot=True, fpass=[]):
    """
    Parameters:
        x (array): Input signal.
        fs (float): Sampling frequency.
        window_length (int): Length of the window.
        OverlapLength (int): Overlap length between windows.
        FFTL (int): FFT length.
        plot (bool): Whether to plot the STFT spectrum.
        fpass (list): Frequency range [fmin, fmax] to display.
    
    Returns:
        s (array): STFT result (complex values).
        f (array): Frequency bins.
        t (array): Time bins.
    """
    # Compute STFT
    f, t, s = stft(x, fs, window=hann(window_length), nperseg=window_length, noverlap=OverlapLength, nfft=FFTL)
    
    # Filter frequency range if fpass is provided
    if len(fpass) == 0:
        fpass = [0, fs / 2]
    
    f_start = np.argmax(f >= fpass[0] - 1e-6)  # Find the start index
    f_end = np.argmax(f >= fpass[1] - 1e-6)    # Find the end index
    if f_end == 0:  # If no frequency exceeds fpass[1], use the last index
        f_end = len(f) - 1
    
    f = f[f_start:f_end + 1]
    s = s[f_start:f_end + 1, :]
    
    # Plot the STFT spectrum
    if plot:
        plt.figure(figsize=(10, 6))
        sdb = 20 * np.log10(np.abs(s))  # Convert magnitude to dB
        plt.imshow(sdb[::-1, :],  # Flip frequency axis for correct orientation
                   extent=[t[0], t[-1], f[0], f[-1]],
                   aspect='auto', cmap='turbo',
                   norm=Normalize(vmax=np.max(sdb), vmin=np.max(sdb) - 60))
        plt.colorbar(label='Magnitude (dB)')
        plt.xlabel('Time (s)')
        plt.ylabel('Frequency (Hz)')
        plt.title('STFT Spectrum (dB)')
        plt.gca().set_ylim(fpass)  # Set frequency limits
        plt.gca().tick_params(axis='both', labelsize=12)
        plt.show()
    
    return s, f, t

def MFTFCCA_compute_calibrations(X, fs, downsample_ratio, win, time_start, time_end, maxlag1, maxlag2, real_STFT, FFT_num=None, fpass=None, plot=True):
    """
    Python implementation of the MATLAB function MFTFCCA_compute_calibrations.
    """
    if FFT_num is None:
        FFT_num = win
    
    # Calculate lag limits
    lag_max = int(np.ceil(maxlag1 * fs)) + int(np.ceil((win - 1) / 2))
    lag_max2 = int(np.floor(maxlag2 * fs)) - int(np.ceil((win - 1) / 2))
    
    # Calculate indices for downsampled and original signals
    x22_r1 = int(np.ceil(time_start * fs / downsample_ratio)) + 1
    x22_r2 = int(np.floor(time_end * fs / downsample_ratio)) + 1
    LH_index = lambda x: (x - 1) * downsample_ratio + 1
    s1_r1 = LH_index(x22_r1)
    s1_r2 = LH_index(x22_r2)
    
    # Validate indices
    x_len = X.shape[1]
    assert s1_r1 > 0, "time_start is too small"
    assert s1_r2 <= x_len, "time_end is too large"
    
    if s1_r1 < lag_max:
        lag_max = s1_r1 - 1
        print(f"maxlag1 is too large, set to {lag_max}")
    
    sx_len = x_len - win + 1
    if sx_len - s1_r2 < lag_max2:
        lag_max2 = sx_len - s1_r2
        print(f"maxlag2 is too large, set to {lag_max2}")
    
    # Calculate lags
    lags = (np.arange(-lag_max + 1, lag_max2 + 1) / fs) + (win - 1) / (2 * fs)
    r = np.zeros((X.shape[0], lag_max + lag_max2))
    
    # Downsample the signal
    X_downsampled = X[:, ::downsample_ratio]
    
    # Compute canonical correlations
    for i in range(X.shape[0]):
        x_high = X[i, :]
        x_low = X_downsampled[i, :]
        
        # Compute STFT
        f, t, s1 = STFT(x_high, fs, win, win // 2, FFT_num, plot=False, fpass=fpass)
        s1 = (s1 - np.mean(s1, axis=1, keepdims=True)) / np.std(s1, axis=1, keepdims=True)  # Z-score normalization
        if real_STFT:
            s1 = np.abs(s1)
        
        # Compute canonical correlations for each lag
        for lag in range(-lag_max + 1, lag_max2 + 1):
            s1_lagged = s1[:, s1_r1 + lag - 1::downsample_ratio][:, :x22_r2 - x22_r1 + 1].T
            x_low_segment = x_low[x22_r1 - 1:x22_r2].reshape(-1, 1)
            r[i, lag_max + lag - 1] = canoncorr(s1_lagged, x_low_segment)[2]
    
    # Normalize correlations
    r_cal_max = np.max(r, axis=1, keepdims=True)
    r_cal_min = np.min(r, axis=1, keepdims=True)
    calibrations = np.hstack((r_cal_max, r_cal_min))
    r = (r - r_cal_min) / (r_cal_max - r_cal_min)
    
    # Plot results
    if plot:
        plt.figure(figsize=(8, 6))
        plt.rcParams.update({'font.size': 18})
        if X.shape[0] > 1:
            r_mean = np.mean(r, axis=0)
            r_std = np.std(r, axis=0)
            plt.fill_between(lags, r_mean - r_std, r_mean + r_std, color='red', alpha=0.3, edgecolor='none')
        else:
            r_mean = r[0, :]
        
        plt.plot(lags, r_mean, linewidth=4, color='blue')
        plt.axvline(-win / (2 * fs), color='black', linestyle='--', linewidth=3)
        plt.axvline(win / (2 * fs), color='black', linestyle='--', linewidth=3)
        plt.xlabel("Time lag (s)", fontsize=20)
        plt.ylabel("Can. Corr.", fontsize=20)
        plt.grid(True)
        plt.show()
    
    return calibrations

def MFTFCCA_lagCC(X, Y, fs, downsample_ratio, win, time_start, time_end, maxlag1, maxlag2, 
                  real_STFT, 
                  calibrations=np.array([]), 
                  FFT_num=np.nan, 
                  surrogate_num=np.nan,
                  fpass=[],
                  partial_high_set=[],
                  partial_low_set=[],
                  plot=True):

    ntrials = X.shape[0]
    x_len = X.shape[1]
    if FFT_num is np.nan:
        FFT_num = win
    if surrogate_num is np.nan:
        surrogate_num = np.round(2000 / ntrials)

    if calibrations.size == 0:
        r_cal_max = np.ones(ntrials)
        r_cal_min = np.zeros(ntrials)
    else:
        r_cal_max = calibrations[:, 0]
        r_cal_min = calibrations[:, 1]

    partial_high = len(partial_high_set) > 0
    partial_low = len(partial_low_set) > 0
    if partial_low:
        fs2 = fs / downsample_ratio
        win2 = np.ceil(win / downsample_ratio)
        FFT_num2 = FFT_num / downsample_ratio

    LH_index = lambda x: (x-1) * downsample_ratio + 1
    HL_index = lambda x: np.ceil((x-1) / downsample_ratio) + 1

    x22_r1 = np.ceil(time_start * fs / downsample_ratio) + 1
    x22_r2 = np.floor(time_end * fs / downsample_ratio) + 1
    s1_r1 = LH_index(x22_r1)
    s1_r2 = LH_index(x22_r2)

    lag_max = np.ceil(maxlag1 * fs) + np.ceil((win - 1) / 2)
    lag_max2 = np.floor(maxlag2 * fs) - np.ceil((win - 1) / 2)

    if s1_r1 < lag_max:
        lag_max = s1_r1 - 1
        print("maxlag1 is too large, set to", lag_max)
    sx_len = x_len - win + 1
    if sx_len - s1_r2 < lag_max2:
        lag_max2 = sx_len - s1_r2
        print("maxlag2 is too large, set to", lag_max2)

    lags = (np.arange(-lag_max + 1, lag_max2 + 1) / fs + (win - 1) / 2 / fs)


    # -------- Compute the original data --------

    r = np.zeros((ntrials, lag_max + lag_max2))

    for i in range(ntrials):
        x_high = X[i, :]
        x_low = Y[i, :]

        s1, _, _ = STFT(x_high, fs, win, win - 1, FFT_num, plot=0, fpass=fpass)
        s1 = zscore(s1.T).T
        if real_STFT:
            s1 = np.abs(s1)

        if partial_high:
            s3_high = []
            if isinstance(partial_high_set, list):
                for k in range(len(partial_high_set)):
                    x_partial = partial_high_set[k][i, :]
                    s3_, _, _ = STFT(x_partial, fs, win, win - 1, FFT_num, plot=0, fpass=fpass)
                    s3_ = zscore(s3_.T).T
                    if real_STFT:
                        s3_ = np.abs(s3_)
                    s3_high.append(s3_)
                s3_high = np.concatenate(s3_high, axis=1)
            else:
                x_partial = partial_high_set[i, :]
                s3_, _, _ = STFT(x_partial, fs, win, win - 1, FFT_num, plot=0, fpass=fpass)
                s3_ = zscore(s3_.T).T
                if real_STFT:
                    s3_ = np.abs(s3_)
                s3_high = s3_

        elif partial_low:
            s3_low = []
            if isinstance(partial_low_set, list):
                for k in range(len(partial_low_set)):
                    x_partial = partial_low_set[k]
                    x_partial = x_partial[i, :]
                    s3_, _, _ = STFT(x_partial, fs2, win2, win2-1, FFT_num2, plot=0, fpass=fpass)
                    s3_ = zscore(s3_.T).T
                    if real_STFT:
                        s3_ = np.abs(s3_)
                    s3_low.append(s3_)
                s3_low = np.concatenate(s3_low, axis=1)
            else:
                x_partial = partial_low_set[i, :]
                s3_, _, _ = STFT(x_partial, fs2, win2, win2-1, FFT_num2, plot=0, fpass=fpass)
                s3_ = zscore(s3_.T).T
                if real_STFT:
                    s3_ = np.abs(s3_)
                s3_low = s3_

        for lag in range(-lag_max + 1, lag_max2 + 1):
            s3 = []
            if partial_high:
                s3 = np.concatenate((s3, s3_high[:, s1_r1 + lag: downsample_ratio: s1_r2 + lag]), axis=1)
            if partial_low:
                s3 = np.concatenate((s3, s3_low[:, HL_index(s1_r1 + lag): HL_index(s1_r2 + lag)]), axis=1)

            if partial_high or partial_low:
                _, _, ri = pCCA(s1[:, s1_r1 + lag: downsample_ratio: s1_r2 + lag].T, x_low[x22_r1:x22_r2].T, s3.T)
            else:
                _, _, ri = canoncorr(s1[:, s1_r1 + lag: downsample_ratio: s1_r2 + lag].T, x_low[x22_r1:x22_r2].T)
            r[i, lag_max + lag] = ri

    r = (r - r_cal_min[:ntrials]) / (r_cal_max[:ntrials] - r_cal_min[:ntrials])
    if ntrials > 1:
        r = np.vstack((np.mean(r, axis=0), np.sqrt(np.var(r, axis=0))))  # / np.sqrt(ntrials)

    r_0 = r

    # -------- Compute the surrogate data --------

    r = np.zeros((ntrials, lag_max + lag_max2, surrogate_num))

    for i in range(ntrials):
        x_high = X[i, :]
        x_low = Y[i, :]
        x_high_su = phaseran(x_high.T, surrogate_num)

        if partial_high:
            s3_high = []
            if isinstance(partial_high_set, list):
                for k in range(len(partial_high_set)):
                    x_partial = partial_high_set[k]
                    x_partial = x_partial[i, :]
                    s3_, _, _ = STFT(x_partial, fs, win, win - 1, FFT_num, plot=0, fpass=fpass)
                    s3_ = zscore(s3_.T).T
                    if real_STFT:
                        s3_ = np.abs(s3_)
                    s3_high.append(s3_)
                s3_high = np.concatenate(s3_high, axis=1)
            else:
                x_partial = partial_high_set[i, :]
                s3_, _, _ = STFT(x_partial, fs, win, win - 1, FFT_num, plot=0, fpass=fpass)
                s3_ = zscore(s3_.T).T
                if real_STFT:
                    s3_ = np.abs(s3_)
                s3_high = s3_

        elif partial_low:
            s3_low = []
            if isinstance(partial_low_set, list):
                for k in range(len(partial_low_set)):
                    x_partial = partial_low_set[k]
                    x_partial = x_partial[i, :]
                    s3_, _, _ = STFT(x_partial, fs2, win2, win2-1, FFT_num2, plot=0, fpass=fpass)
                    s3_ = zscore(s3_.T).T
                    if real_STFT:
                        s3_ = np.abs(s3_)
                    s3_low.append(s3_)
                s3_low = np.concatenate(s3_low, axis=1)
            else:
                x_partial = partial_low_set[i, :]
                s3_, _, _ = STFT(x_partial, fs2, win2, win2-1, FFT_num2, plot=0, fpass=fpass)
                s3_ = zscore(s3_.T).T
                if real_STFT:
                    s3_ = np.abs(s3_)
                s3_low = s3_

        for s in range(surrogate_num):
            x_high = x_high_su[:, :, s].T
            x_low = phaseran(x_low.T, 1).T

            s1, _, _ = STFT(x_high, fs, win, win - 1, FFT_num, plot=0, fpass=fpass)

            s1 = zscore(s1.T).T
            if real_STFT:
                s1 = np.abs(s1)

            for lag in range(-lag_max + 1, lag_max2 + 1):
                s3 = []
                if partial_high:
                    s3 = np.concatenate((s3, s3_high[:, s1_r1 + lag: downsample_ratio: s1_r2 + lag]), axis=1)
                if partial_low:
                    s3 = np.concatenate((s3, s3_low[:, HL_index(s1_r1 + lag): HL_index(s1_r2 + lag)]), axis=1)

                if partial_high or partial_low:
                    _, _, ri = pCCA(s1[:, s1_r1 + lag: downsample_ratio: s1_r2 + lag].T, x_low[x22_r1:x22_r2].T, s3.T)
                else:
                    _, _, ri = canoncorr(s1[:, s1_r1 + lag: downsample_ratio: s1_r2 + lag].T, x_low[x22_r1:x22_r2].T)

                r[i, lag_max + lag, s] = ri

    r = (r - r_cal_min) / (r_cal_max - r_cal_min)
    r = np.vstack((np.mean(r, axis=(0, 2)), np.sqrt(np.var(r, axis=(0, 2)))))  # / np.sqrt(ntrials)
    r_su = r


    if plot:

        fig = plt.figure(figsize=(8, 4))
        ax = fig.add_subplot(111)

        norm_factor = np.sqrt(ntrials)

        r = r_su
        ax.fill_between(lags, r[0, :] + r[1, :] / norm_factor, r[0, ::-1] - r[1, ::-1] / norm_factor,
                        facecolor='b', edgecolor='none')

        if ntrials > 1:
            r = r_0
            ax.fill_between(lags, r[0, :] + r[1, :] / norm_factor, r[0, ::-1] - r[1, ::-1] / norm_factor,
                            facecolor='r', edgecolor='none')

        ax.plot(lags, r_su[0, :], linewidth=6, color='b')
        ax.plot(lags, r_0[0, :], linewidth=6, color='r')

        y_limits = ax.get_ylim()

        ax.plot([-win / fs / 2, -win / fs / 2], y_limits, 'k--', linewidth=3)
        ax.plot([win / fs / 2, win / fs / 2], y_limits, 'k--', linewidth=3)

        ax.legend(["Surrogate", "Original"], loc="best")
        ax.set_xlabel("Time lag (s)")
        ax.set_ylabel("Can. Corr")
        ax.tick_params(labelsize=25)
        ax.spines['linewidth'] = 4

        plt.show()

    return lags, r_0, r_su

def MFTFCCA_spectral(X, Y, fs, downsample_ratio, 
                     win, time_start, time_end,
                     real_STFT,
                     calibrations=np.array([]),
                     lag=0,
                     FFT_num=np.nan,
                     surrogate_num=np.nan,
                     fpass=[],
                     partial_high_set=[],
                     partial_low_set=[],
                     plot=True,
                     plot_interval=np.nan,):
    
    ntrials, x_len = X.shape
    if FFT_num is np.nan:
        FFT_num = win
    if surrogate_num is np.nan:
        surrogate_num = np.round(2000 / ntrials)
    if plot_interval is np.nan:
        plot_interval = np.round(ntrials / 100)
    
    if calibrations.size == 0:
        r_cal_max = np.ones(ntrials)
        r_cal_min = np.zeros(ntrials)
    else:
        r_cal_max = calibrations[:, 0]
        r_cal_min = calibrations[:, 1]

    partial = len(partial_high_set) > 0
    partial_low = len(partial_low_set) > 0
    if partial_low:
        fs2 = fs / downsample_ratio
        win2 = np.ceil(win / downsample_ratio)
        FFT_num2 = FFT_num / downsample_ratio

    LH_index = lambda x: (x-1) * downsample_ratio + 1

    x22_r1 = np.ceil(time_start * fs / downsample_ratio) + 1
    x22_r2 = np.floor(time_end * fs / downsample_ratio) + 1
    lag = round(lag * fs - (win - 1) / 2)
    s1_r1 = LH_index(x22_r1) + lag
    s1_r2 = LH_index(x22_r2) + lag

    assert s1_r1 >= 1, "time_start is too small for the lag"
    assert x_len - win + 1 >= s1_r2, "time_end is too large for the lag"

    As0 = []
    r_ori = np.zeros(ntrials)

    for i in range(ntrials):
        x_high = X[i, :]
        x_low = Y[i, :]

        _, _, s1 = stft(x_high, fs, window='hann', nperseg=win, noverlap=win-1, nfft=FFT_num, boundary=None)
        s1 = (s1.T - np.mean(s1, axis=1)).T
        if real_STFT:
            s1 = np.abs(s1)
        s1 = s1[:, int(s1_r1):int(s1_r2):int(downsample_ratio)]

        if partial:
            x_partial = partial_high_set[i, :]
            _, _, s3 = stft(x_partial, fs, window='hann', nperseg=win, noverlap=win-1, nfft=FFT_num, boundary=None)
            s3 = (s3.T - np.mean(s3, axis=1)).T
            if real_STFT:
                s3 = np.abs(s3)
            s3 = s3[:, int(s1_r1):int(s1_r2):int(downsample_ratio)]
            A, _, ri = pCCA(s1.T, x_low[int(x22_r1):int(x22_r2)].T, s3.T)
        elif partial_low:
            x_partial = partial_low_set[i, :]
            _, _, s3 = stft(x_partial, fs2, window='hann', nperseg=win2, noverlap=win2-1, nfft=FFT_num2, boundary=None)
            s3 = (s3.T - np.mean(s3, axis=1)).T
            if real_STFT:
                s3 = np.abs(s3)
            s3 = s3[:, int(HL_index(s1_r1)):int(HL_index(s1_r2)):int(downsample_ratio)]
            A, _, ri = pCCA(s1.T, x_low[int(x22_r1):int(x22_r2)].T, s3.T)
        else:
            _, _, ri = CCA().fit(s1.T, x_low[int(x22_r1):int(x22_r2)].T).transform(s1.T, x_low[int(x22_r1):int(x22_r2)].T)

        ri = (ri - r_cal_min[i]) / (r_cal_max[i] - r_cal_min[i])
        A = A / np.mean(np.abs(A))

        As0.append(A)
        r_ori[i] = ri

    if plot:
        plot_CCA_coeff(stft_f, As0, plot_interval=plot_interval)

    As = []
    r = np.zeros((ntrials, 2))
    r_shuffle = np.zeros((ntrials, surrogate_num))
    for i in range(ntrials):
        x_high = X[i, :]
        x_low = Y[i, :]
        x_high_set = phaseran(x_high.T, surrogate_num)
        r_j = np.zeros(surrogate_num)

        for j in range(surrogate_num):
            x_high = x_high_set[:, :, j].T
            x_low = phaseran(x_low.T, 1).T
            _, _, s1 = stft(x_high, fs, window='hann', nperseg=win, noverlap=win-1, nfft=FFT_num, boundary=None)
            s1 = (s1.T - np.mean(s1, axis=1)).T
            if real_STFT:
                s1 = np.abs(s1)
            s1 = s1[:, int(s1_r1):int(s1_r2):int(downsample_ratio)]

            if partial:
                x_partial = partial_high_set[i, :]
                _, _, s3 = stft(x_partial, fs, window='hann', nperseg=win, noverlap=win-downsample_ratio, nfft=FFT_num, boundary=None)
                s3 = (s3.T - np.mean(s3, axis=1)).T
                if real_STFT:
                    s3 = np.abs(s3)
                s3 = s3[:, int(s1_r1):int(s1_r2):int(downsample_ratio)]
                A, _, ri = pCCA(s1.T, x_low[int(x22_r1):int(x22_r2)].T, s3.T)
            elif partial_low:
                x_partial = partial_low_set[i, :]
                _, _, s3 = stft(x_partial, fs2, window='hann', nperseg=win2, noverlap=win2-1, nfft=FFT_num2, boundary=None)
                s3 = (s3.T - np.mean(s3, axis=1)).T
                if real_STFT:
                    s3 = np.abs(s3)
                s3 = s3[:, int(HL_index(s1_r1)):int(HL_index(s1_r2)):int(downsample_ratio)]
                A, _, ri = pCCA(s1.T, x_low[int(x22_r1):int(x22_r2)].T, s3.T)
            else:
                _, _, ri = CCA().fit(s1.T, x_low[int(x22_r1):int(x22_r2)].T).transform(s1.T, x_low[int(x22_r1):int(x22_r2)].T)

            ri = (ri - r_cal_min[i]) / (r_cal_max[i] - r_cal_min[i])
            A = A / np.mean(np.abs(A))
            As.append(A)

            r_j[j] = ri

        r[i, 0] = np.mean(r_j)
        r[i, 1] = np.var(r_j)

        r_shuffle[i, :] = r_j

    t_stat, p_value = ttest_ind(r_ori, r[:, 0], equal_var=False)
    if t_stat > 0:
        print(f"ttest: Causal, p={p_value:.2e}")
    else:
        print(f"ttest: Non-causal, p={p_value:.2e}")

    ks_stat, p_value = ks_2samp(r_ori, r[:, 0])
    if ks_stat > 0:
        print(f"KStest: Causal, p={p_value:.2e}")
    else:
        print(f"KStest: Non-causal, p={p_value:.2e}")

    return stft_f, r_ori, As0, r_shuffle, As