%% G_cause
% 
% Estimate the Granger causality (time / frequency domain)
% 
%% Arguments
%
% _input_
%
%     X          variate matrix (nvar * nobs * ntrials)
%     fs         sampling rate (Hz)
%     options
%       seed     random seed (default: 0, unseed)
%       plot     plot GC (default: false)
%       mode     GC estimation mode. {t, f, tf} (default: tf)
%                t: in time domain  f: in freq domain  tf: both domains
%       report   report all estimated params (default: true)
% 
% _output_
%
%     GC_f       spectral conditional GC (nvar * nvar * nf)
%     f          frequencies corresponding to GC_f (1 * nf)
%
%%
function [GC_f, f] = G_cause(X, fs, options)

arguments
    X
    fs
    options.seed (1,1) {mustBeNumeric} = 0
    options.plot (1,1) = false
    options.mode (1,1) = "tf"
    options.report (1,1) = true
    options.maxorder (1,1) = 30

end

seed = options.seed;
plot = options.plot;
mode = options.mode;
report = options.report;
momax = options.maxorder;

rng_seed(seed);

regmode   = 'LWR';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
% momax     = 50;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = 'F';    % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDRD'; % multiple hypothesis test correction (see routine 'significance')

fres      = [1024];    % frequency resolution (empty for automatic calculation)

[AIC2, BIC2, moAIC2, moBIC2] = tsdata_to_infocrit(X, momax, icregmode, report);

morder = moAIC2;
[A, SIG] = tsdata_to_var(X, morder, regmode);
assert(~isbad(A),'VAR estimation failed - bailing out');
info = var_info(A, SIG, report);

if report
    fprintf('A:')
    A
    fprintf('SIG:')
    SIG
end

% time
if mode == "t" || mode == "tf"
    [F,pval] = var_to_pwcgc(A,SIG,X,regmode,tstat);
    assert(~isbad(F,false),'GC calculation failed - bailing out');
    sig = significance(pval,alpha,mhtc);
    if plot
        figure(2); clf;
        sgtitlex('Pairwise-conditional Granger causality - time domain');
        subplot(1,3,1);
        plot_pw(F);
        title('Pairwise-conditional GC');
        subplot(1,3,2);
        plot_pw(pval);
        title(['p-values (' tstat '-test)']);
        subplot(1,3,3);
        plot_pw(sig);
        title(['Significant at \alpha = ' num2str(alpha)]);
        hold(gca, 'off')
    end
    if report        
        F
        pval
        sig
    end

end


% freq
if mode == "f" || mode == "tf"
    if isempty(fres)
        fres = 2^nextpow2(info.acdec); % based on autocorrelation decay; alternatively, you could try fres = 2^nextpow2(nobs);
        fprintf('\nfrequency resolution auto-calculated as %d (increments ~ %.2gHz)\n',fres,fs/2/fres);
    end
    if fres > 20000 % adjust to taste
        fprintf(2,'\nWARNING: large frequency resolution = %d - may cause computation time/memory usage problems\nAre you sure you wish to continue [y/n]? ',fres);
        istr = input(' ','s'); if isempty(istr) || ~strcmpi(istr,'y'); fprintf(2,'Aborting...\n'); return; end
    end

    % Calculate spectral pairwise-conditional causalities at given frequency
    % resolution by state-space method.

    ptic('\n*** var_to_spwcgc... ');
    GC_f = var_to_spwcgc(A,SIG,fres);
    f = linspace(0, fs/2, size(GC_f, 3));
    assert(~isbad(GC_f,false),'spectral GC calculation failed - bailing out');
    ptoc;

    % Plot spectral causal graph.

    if plot
        figure();
        sgtitlex('Pairwise-conditional Granger causality - frequency domain');
        plot_spw(GC_f,fs);
    end

end


end