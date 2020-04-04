%% Fig3_heat_markoverr.m
% 
% This code generates Figure 5.3 from the paper
%   'Efficient Randomized Algorithms for
%    Subspace System Identification'
%       -Minster, Saibaba, Kar, Chakrabortty


%load matrices
addpath('test_matrices')

load steel_A.mat
load steel_B.mat
load steel_C.mat
load steel_E.mat


% Transform system into standard form
L = chol(E,'lower');
A = L\(A/(L'));
B = L\B;
C = C/(L');
disp('cholesky')


% parameters
ts = .007;
ni = size(B,2);
no = size(C,1);
p = 1000;
q = p;
N = 2*p-1;
r = 20;


% convert to discrete
sys = ss(A,B,C,0);
sysd = c2d(sys,ts,'tustin');
[Ad,Bd,Cd,Dd] = ssdata(sysd);

%% Markov parameters
markov = cell(1,N);
f = Bd;
for jj = 1:N
    g = Cd*f;
    f = Ad*f;           
    markov{jj} = g;
end

%% System ID Algorithms

% Full SVD-ERA
[Arf,Brf,Crf,~] = impulse_era(markov,s,l,m,r,'full');

% SVDS-H
[Ar_svds,Br_svds,Cr_svds,~] = impulse_era(markov,s,l,m,r,'svdshankel'); 

% RandSVD-H
[Ar_rsvdh,Br_rsvdh,Cr_rsvdh,~] = impulse_era(markov,s,l,m,r,'randsvdhankel'); 


%% Form new markov parameters
markov_f = cell(1,N);
f = Brf;
for jj = 1:N
    g = Crf*f;
    f = Arf*f;           
    markov_f{jj} = g;
end

markov_svds = cell(1,N);
f = Br_svds;
for jj = 1:N
    g = Cr_svds*f;
    f = Ar_svds*f;           
    markov_svds{jj} = g;
end

markov_rsvdh = cell(1,N);
f = Br_rsvdh;
for jj = 1:N
    g = Cr_rsvdh*f;
    f = Ar_rsvdh*f;           
    markov_rsvdh{jj} = g;
end

%% Compute relative error in Markov parameters
err_svds = zeros(1,N);
err_rsvdh = zeros(1,N);
for k = 1:N
    err_svds(k) = norm(markov_f{k}-markov_svds{k})/norm(markov_f{k});
    err_rsvdh(k) = norm(markov_f{k}-markov_rsvdh{k})/norm(markov_f{k});
end

%% plot error
figure,
semilogy(1:N,err_svds,'linewidth',2), hold on
semilogy(1:N,err_rsvdh,'linewidth',2)
legend('SVDS-H','RandSVD-H')
title('Markov Parameter Error')
ylabel('Relative Error')
xlabel('$k$','interpreter','latex')
set(gca,'fontsize',18)



