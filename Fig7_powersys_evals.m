%% Fig7_powersys_evals.m
% 
% This code generates Figure 5.7 from the paper
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

%% System ID Algorithms

% RandSVD-H
[Ar_rsvdh,Br_rsvdh,Cr_rsvdh,~] = impulse_era(markov,s,l,m,r,'randsvdhankel'); 

% RandTERA
[Ar_rtera,Br_rtera,Cr_rtera,~] = impulse_era(markov,s,l,m,r,'randtera');

%% Plot eigenvalues

figure,
subplot(1,3,1)
plot(real(eig(Ad)),imag(eig(Ad)),'bo','markersize',10)
xlim([.6,1])
title('Original')
set(gca,'fontsize',18)

subplot(1,3,2)
plot(real(eig(Ar_rsvdh)),imag(eig(Ar_rsvdh)),'k*','markersize',10)
title('RandSVD-H')
set(gca,'fontsize',18)

subplot(1,3,3)
plot(real(eig(Ar_rtera)),imag(eig(Ar_rtera)),'rx','markersize',10)
title('RandTERA')
set(gca,'fontsize',18)