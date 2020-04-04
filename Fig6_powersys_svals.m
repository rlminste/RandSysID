%% Fig6_powersys_svals.m
% 
% This code generates Figure 5.6 from the paper
%   'Efficient Randomized Algorithms for
%    Subspace System Identification'
%       -Minster, Saibaba, Kar, Chakrabortty


%load matrices
addpath('test_matrices')
load powersystem155.mat

%parameters
ts = .007;
m = size(B,2);
l = size(C,1);
s = 300;
N = 2*p-1;
r = size(A,1); %no model reduction

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

%% computing singular values

% SVD
H = hankelize(markov,s,l,m);
sf = svd(H,'econ');

% RandSVD-H
rho = 20;
maxiter = 1;

Hf = makehankelfun(markov,s,l,m);
[~,S_rsvdh,~] = rsvdFun(Hf,s*m,r,rho,maxiter);
s_rsvdh = diag(S_rsvdh);

% TERA
epsilon = [.1,.05,.01];
s_tera = cell(1,3);
for i = 1:3
    [markov_proj, W1, W2,~] = tera(markov,epsilon(i));     
    l1 = size(W1,2); 
    m1 = size(W2,2); 
    H = hankelize(markov_proj,s,l1,m1);
    [~,S,~] = svd(H,0);
    s_tera{i} = diag(S);
end
s1_tera = s_tera{1};
s2_tera = s_tera{2};
s3_tera = s_tera{3};


% RandTERA
s_rtera = cell(1,3);
for i = 1:3
    [markov_proj, W1, W2,~] = tera(markov,epsilon(i));  %W1 is no x r1, W2 is ni x r2 
    l1 = size(W1,2);
    m1 = size(W2,2);
    H = markov_proj(1:N);
    Hf = makehankelfun(H,s,l1,m1,0);
    [~,S,~] = rsvdFun(Hf,r,rho,maxiter);

    s_rtera{i} = diag(S);
end


%% plot singular values
figure,

%singular values with TERA
subplot(1,2,1)
semilogy(1:155,sf(1:155),'k','linewidth',2), hold on
semilogy(1:155,s_rsvdh,'r--','linewidth',2)
semilogy(1:155,s1_tera(1:155),'b-.','linewidth',2)
semilogy(1:155,s2_tera(1:155),'m--','linewidth',2.5)
semilogy(1:155,s3_tera(1:155),'g:','linewidth',2.5)
legend('Full SVD','RandSVD-H','TERA: .1','TERA: .05','TERA: .01','location','best')
title('Singular Values of H')
set(gca,'fontsize',18)

%singular values with RandTERA
subplot(1,2,2)
semilogy(1:155,sf(1:155),'k','linewidth',2), hold on
semilogy(1:length(srsvdh),srsvdh,'r--','linewidth',2)
semilogy(1:length(s_rtera{1}),s_rtera{1},'b-.','linewidth',2)
semilogy(1:length(s_rtera{2}),s_rtera{2},'m:','linewidth',2.5)
semilogy(1:length(s_rtera{3}),s_rtera{3},'g:','linewidth',2.5)
legend('Full SVD','RandSVD-H','RandTERA: .1','RandTERA: .05','RandTERA: .01','location','best')
title('Singular Values of H')
set(gca,'fontsize',18)


