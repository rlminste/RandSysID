function [markov_proj, W1, W2,t_newmarkov] = tera(markov,epsilon)
% 
% Implements steps 1-7 of algorithm 1 (Gugercin and Kramer)
% 
% Inputs:
%    - markov: cell (1 x N) of markov paramers h_i (l x m) i = 1,...,N
%    - epsilon: threshold to determine how many singular values to keep
% 
% Outputs:
%    - markov_proj: cell (1 x N) of markov paramers h_i (r1 x r2) i = 1,...,N
%    - W1: projection basis l x r1 
%    - W2: projection basis m x r2
%    - t_newmarkov: time of forming new markov parameters
%
% - Written by Arvind K. Saibaba and Rachel Minster, 2020
%
% Reference: Kramer and Gugercin, 'Tangential and Interpolation-based
%       Eigensystem Realization Algorithm for MIMO Systems', Mathematical and
%       Computer Modelling of Dynamical Systems, 2016


N = length(markov);

% Compute W1
H = cell2mat(markov);
[U1,S1,~] = svd(H,'econ');

s1 = diag(S1);
s1 = s1./s1(1);
ind = find(s1<epsilon,1);
r1 = ind-1;
   
    
% Compute W2 
markovt = cell(N,1);
for j = 1:N
    markovt{j} = markov{j};
end

Ht = cell2mat(markovt);
[~,S2,V2] = svd(Ht,0);

s2 = diag(S2);
s2 = s2./s2(1);
ind = find(s2<epsilon,1);
r2 = ind-1;
    
% Compute bases used to compute projections
W1 = U1(:,1:r1);
W2 = V2(:,1:r2);

% Compute new Markov parameters
tic; markov_proj = cell(1,N);
for j = 1:N
   markov_proj{j} = W1'*(markov{j}*W2); 
end
t_newmarkov = toc;

end