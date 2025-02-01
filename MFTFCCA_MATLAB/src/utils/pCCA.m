function[ccaEigvector1, ccaEigvector2, D] = pCCA(X, Y, Z)


if nargin == 2

    X = X - mean(X,1);
    Y = Y - mean(Y,1);
    
    dataLen1 = size(X, 2);
    dataLen2 = size(Y, 2);
    data = [X Y];

    covariance = cov(data);

    % Sxx = covariance(1 : dataLen1, 1 : dataLen1) + eye(dataLen1) * 10^(-7);

    Sxx = covariance(1 : dataLen1, 1 : dataLen1);

    % Syy = covariance(dataLen1 + 1 : size(covariance, 2), dataLen1 + 1 : size(covariance, 2)) ...

    % + eye(dataLen2) * 10^(-7);

    Syy = covariance(dataLen1 + 1 : dataLen1 + dataLen2, dataLen1 + 1 : dataLen1 + dataLen2);

    Sxy = covariance(1 : dataLen1, dataLen1 + 1 : dataLen1 + dataLen2);

    % using SVD to compute the projection

    Hx = (Sxx)^(-1/2);

    Hy = (Syy)^(-1/2);


    H = Hx * Sxy * Hy;

    [U, D, V] = svd(H, 'econ');

    ccaEigvector1 = Hx * U;

    ccaEigvector2 = Hy * V;

    % make the canonical correlation variable has unit variance

    ccaEigvector1 = ccaEigvector1 * diag(diag((eye(size(ccaEigvector1, 2)) ./ sqrt(ccaEigvector1' * Sxx * ccaEigvector1))));

    ccaEigvector2 = ccaEigvector2 * diag(diag((eye(size(ccaEigvector2, 2)) ./ sqrt(ccaEigvector2' * Syy * ccaEigvector2))));

    
else

    % Center the variables
    X = X - mean(X,1);
    Y = Y - mean(Y,1);
    Z = Z - mean(Z,1);

    dataLen1 = size(X, 2);

    dataLen2 = size(Y, 2);

    dataLen3 = size(Z, 2);


    % Construct the scatter of each view and the scatter between them

    data = [X Y Z];

    covariance = cov(data);

    % Sxx = covariance(1 : dataLen1, 1 : dataLen1) + eye(dataLen1) * 10^(-7);

    Sxx = covariance(1 : dataLen1, 1 : dataLen1);

    % Syy = covariance(dataLen1 + 1 : size(covariance, 2), dataLen1 + 1 : size(covariance, 2)) ...

    % + eye(dataLen2) * 10^(-7);

    Syy = covariance(dataLen1 + 1 : dataLen1 + dataLen2, dataLen1 + 1 : dataLen1 + dataLen2);

    Sxy = covariance(1 : dataLen1, dataLen1 + 1 : dataLen1 + dataLen2);

    % Syx = Sxy';

    Szz = covariance(dataLen1+dataLen2+1:end, dataLen1+dataLen2+1:end);

    Sxz = covariance(1:dataLen1, dataLen1+dataLen2+1:end);
    Szx = Sxz';
    Syz = covariance(dataLen1 + 1 : dataLen1 + dataLen2, dataLen1+dataLen2+1:end);
    Szy = Syz';

    % Sxx = cov(data1, data1);
    % Syy = cov(data2, data2);
    % Szz = cov(data3, data3);
    % Sxz = cov(data1, data3);
    % Szx = Sxz';
    % Syz = cov(data2, data3);
    % Szy = Syz';

    % Szz1 = (Szz+eye(length(Szz))*1e-5)^(-1);
    Szz1 = Szz ^ (-1);

    Sxxz = Sxx - Sxz * Szz1 * Szx;
    Syyz = Syy - Syz * Szz1 * Szy;
    Sxyz = Sxy - Sxz * Szz1 * Szy;

    Sxx = Sxxz;
    Syy = Syyz;
    Sxy = Sxyz;


    % using SVD to compute the projection

    Hx = (Sxx)^(-1/2);

    Hy = (Syy)^(-1/2);


    H = Hx * Sxy * Hy;

    [U, D, V] = svd(H, 'econ');

    ccaEigvector1 = Hx * U;

    ccaEigvector2 = Hy * V;

    % make the canonical correlation variable has unit variance

    ccaEigvector1 = ccaEigvector1 * diag(diag((eye(size(ccaEigvector1, 2)) ./ sqrt(ccaEigvector1' * Sxx * ccaEigvector1))));

    ccaEigvector2 = ccaEigvector2 * diag(diag((eye(size(ccaEigvector2, 2)) ./ sqrt(ccaEigvector2' * Syy * ccaEigvector2))));

end

end