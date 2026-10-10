function [ccaEigvector1, ccaEigvector2, D] = pCCA(X, Y, Z, ridge)
% pCCA  (Partial) canonical correlation analysis.
%
%   [A, B, D] = pCCA(X, Y)        canonical correlation between X and Y
%   [A, B, D] = pCCA(X, Y, Z)     canonical correlation between X and Y
%                                 given Z (partial CCA)
%   [A, B, D] = pCCA(X, Y, Z, ridge)   as above, with relative Tikhonov
%                                 regularisation of the residual
%                                 covariances (default 0).
%
%   Rows are observations, columns are variables. Z may contain any number
%   of columns; they are conditioned on jointly. Redundant (collinear or
%   duplicated) columns of Z are handled without error and contribute
%   nothing, because X and Y are projected onto the orthogonal complement
%   of span(Z) rather than through an explicit inv(Szz).

if nargin < 4 || isempty(ridge)
    ridge = 0;
end

X = X - mean(X,1);
Y = Y - mean(Y,1);

% Variance scale before partialling out, used as the reference for the
% rank tolerance below.
scaleX = max(var(X, 0, 1));
scaleY = max(var(Y, 0, 1));

if nargin >= 3 && ~isempty(Z)
    Z = Z - mean(Z,1);

    % Residualise X and Y on the column space of Z. Using the projector
    % Z*pinv(Z) instead of Sxz*inv(Szz)*Szx makes the result invariant to
    % rank deficiency in Z (duplicated / collinear partial-out variables)
    % and avoids squaring the condition number of Z.
    Zp = pinv(Z);
    X = X - Z * (Zp * X);
    Y = Y - Z * (Zp * Y);
end

dataLen1 = size(X, 2);
dataLen2 = size(Y, 2);

covariance = cov([X Y]);

Sxx = covariance(1 : dataLen1, 1 : dataLen1);
Syy = covariance(dataLen1 + 1 : dataLen1 + dataLen2, dataLen1 + 1 : dataLen1 + dataLen2);
Sxy = covariance(1 : dataLen1, dataLen1 + 1 : dataLen1 + dataLen2);

if ridge > 0
    Sxx = Sxx + ridge * trace(Sxx) / dataLen1 * eye(dataLen1);
    Syy = Syy + ridge * trace(Syy) / dataLen2 * eye(dataLen2);
end

% Inverse square roots restricted to the numerically non-null subspace of
% each whitening matrix; directions with no variance left after
% partialling out Z are dropped instead of being inverted.
nobs = size(X, 1);
dimZ = 0;
if nargin >= 3 && ~isempty(Z)
    dimZ = size(Z, 2);
end
Hx = invsqrt_psd(Sxx, scaleX, "X", nobs, dimZ);
Hy = invsqrt_psd(Syy, scaleY, "Y", nobs, dimZ);

H = Hx' * Sxy * Hy;

[U, D, V] = svd(H, 'econ');

ccaEigvector1 = Hx * U;
ccaEigvector2 = Hy * V;

% make the canonical correlation variable has unit variance
ccaEigvector1 = ccaEigvector1 * diag(diag((eye(size(ccaEigvector1, 2)) ./ sqrt(ccaEigvector1' * Sxx * ccaEigvector1))));
ccaEigvector2 = ccaEigvector2 * diag(diag((eye(size(ccaEigvector2, 2)) ./ sqrt(ccaEigvector2' * Syy * ccaEigvector2))));

end


function H = invsqrt_psd(S, scale, side, nobs, dimZ)
% Symmetric inverse square root, truncated at the numerical rank. Variance
% directions that are negligible relative to SCALE (the variance before
% partialling out) are dropped rather than inverted, so that over-
% partialling shows up as lost dimensions instead of amplified noise.
S = (S + S') / 2;
[V, E] = eig(S);
e = real(diag(E));
tol = max(size(S)) * eps(max([e; scale]));
keep = e > max(tol, 0);
if ~any(keep)
    if any(~isfinite(e))
        error("pCCA:nonFinite", ...
              "%s has non-finite values (NaN/Inf) - check the input, not the partialling.", side);
    end
    error("pCCA:rankDeficient", ...
          "No variance left in %s after partialling out Z: %d observations, " + ...
          "dim(%s) = %d, dim(Z) = %d. Z spans the data " + ...
          "(dim(Z) >= observations - 1, or Z contains %s itself).", ...
          side, nobs, side, size(S, 1), dimZ, side);
end
H = V(:, keep) * diag(1 ./ sqrt(e(keep)));
end
