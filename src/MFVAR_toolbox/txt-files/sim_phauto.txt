function X = sim_phauto(E, Bh, B1, p, h)
% PURPOSE: Simulate a (p,h)-autoregression
%--------------------------------------------------------------------------
% USAGE: X = sim_phauto(E, Bh, B1, p, h)
% where: E    = error term (T-by-K matrix)
%        T    = # of observations
%        K    = dimension of the VAR process
%        Bh   = pK-by-K coefficient matrix with predction horizon h
%        B1   = pK-by-K coefficient matrix with predction horizon 1       
%        p    = lag length
%        h    = prediction horizon
%--------------------------------------------------------------------------
% RETURNS: X = simulated process (T-by-K) 
%--------------------------------------------------------------------------
% NOTE: k-th block of Bh contains the transpose of A_k(h) for
%       k = 1,...,p
%--------------------------------------------------------------------------
% References: E.Ghysels, J.B.Hill, & K.Motegi (2013)
%             Testing for Granger Causality with Mixed Frequency Data. 
%             Working Paper at UNC. 
%             Eq. (2.5) gives a (p,h)-autoregression.
% -------------------------------------------------------------------
% Written by Kaiji Motegi, UNC Chapel Hill.
% Jun. 27, 2013.
% -------------------------------------------------------------------------

[T, K] = size(E);



% When h = 1, use Bh as B1.
% In causality_test.m which implements bootstrap, the null of non-causality
% is imposed on Bh only.
if h == 1
    B1 = Bh;
end;    

% Compute impulse response coefficients based on B1 
Psimat = zeros(K, K, h);
Psimat(:,:,1) = eye(K);
if h > 1
    for k = 2:h
        for s = 1:p
              if k > s
                  first = (s-1) * K + 1;
                  last = s * K;
                  Psimat(:,:,k) = Psimat(:,:,k) + B1(first:last,:)' * Psimat(:,:,k-s);
              else
                  break;
              end;
         end;     
    end;    
end;

X = zeros(T+p+h-1, K);
E_ = [zeros(p+h-1, K); E];

% construct X
for t = (p+h):(T+p+h-1)
     % error part
     for k = 1:h
          X(t,:) = X(t,:) + E_(t+1-k,:) * Psimat(:,:,k)';
     end;
     % autoregressive part
     for k = 1:p
          first = (k-1) * K + 1;
          last = k * K;
          X(t,:) = X(t,:) + X(t-h+1-k,:) * Bh(first:last, :);
     end;    
end;     

% cut initial values
X = X((p+h):end, :);


