function [U,S,V] = rsvdFun(Afun,n,k,p,maxiter)
% Computes randomized SVD of a matrix
%
% Inputs: 
%    - Afun: function representing matrix product and transpose
%    - n: size(A,2)
%    - k: target rank
%    - p: oversampling parameter
%    - maxiter: number of subspace iterations
% 
% Outputs: U,S,V so that A = USV'
%
% - Written by Eric Hallman, 2020
% - Modified by Rachel Minster, 2020

Omega = randn(n,k+p);

%% Initial iteration
Y = Afun(Omega,'notransp'); 
[Q,~] = qr(Y,0);

%% Subsequent iterations
    for j = 1:maxiter
        Y = Afun(Q,'transp');
        [Q,~]=qr(Y,0);
        Y = Afun(Q,'notransp');
        [Q,~]=qr(Y,0);
    end
    
%% Low rank approximation
    B = Afun(Q,'transp');
    B = B';
    
    [U,S,V] = svd(B,'econ');
    U = Q*U;
    
 %% Compress
    U = U(:,1:k);
    V = V(:,1:k);
    S = S(1:k,1:k);
end