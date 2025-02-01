function X = sim_VAR(E, Amat)
% PURPOSE: Simulate a VAR(p) process
%--------------------------------------------------------------------------
% USAGE: X = sim_VAR(E, Amat)
% where: E    = error term (T-by-K matrix)
%        T    = # of obs.
%        K    = dimension of the VAR process
%        Amat = K-by-K-by-p coefficient matrix
%               (Amat(:,:,i) = VAR(i) coefficient for i = 1, ... , p)
%        p    = lag length
%--------------------------------------------------------------------------
% RETURNS: X = simulated process (T-by-K) 
%--------------------------------------------------------------------------
% MODEL: X(t,:)' = sum_{i=1}^{p} Amat(:,:,i) * X(t-i,:)' + E(t,:)'
%        X((1-p):0, :) = zeros(p, K);
%--------------------------------------------------------------------------
% Written by Kaiji Motegi, UNC Chapel Hill.
% Aug. 4, 2012.
% -------------------------------------------------------------------------

% T is the number of observations
% K is the dimension of VAR process
[T, K] = size(E);

% lag length
p = size(Amat,3);

% include initial values
X = zeros(T + p, K);

% simulate the VAR(p) process
for t = (p+1):(T+p)
     X(t,:) = E(t-p,:);
     for i = 1:p
          X(t,:) = X(t,:) + X(t-i,:) * Amat(:,:,i)';
     end;
end;    

% cut initial values
X = X((p+1):end,:);
