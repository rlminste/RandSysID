function [Ar,Br,Cr] = era(U,S,V,r,m,l)
%
% Identifies system matrices Ar,Br,Cr
%
% Inputs:
%    - U,S,V: an svd/randomized svd of block Hankel matrix
%    - r: target system size
%    - m: number of inputs
%    - l: number of outputs
%
% Outputs:
%    - Ar,Br,Cr identified discrete-time system matrices
%
% - Originally written by Boris Kramer, 2016
% - Modified by Rachel Minster, 2020

nrH = max(size(U));
ncH = max(size(V));

% Truncate and partition matrices
Uf = U(1:nrH-l, 1:r);   % Cutting last block-row
Ul = U(l+1:nrH, 1:r);   % Cutting first block-row
V1 = V(1:ncH-m , 1:r);
S = S(1:r,1:r);


Uhat = Uf'*Ul;
Ar = Uhat;
B = S*V1';
Br = B(:,1:m);
Cr = Uf(1:l,:).*(ones(l,1));  

end