%% Fig2_heat_evals.m
% 
% This code generates Figure 5.2 from the paper
%   'Efficient Randomized Algorithms for
%    Subspace System Identification'
%       -Minster, Saibaba, Kar, Chakrabortty


%load matrices
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
m = size(B,2);
l = size(C,1);
s = 1000;
N = 2*s-1;
r = 20;


% convert to discrete
sys = ss(A,B,C,0);
sysd = c2d(sys,ts,'tustin');
[Ad,Bd,Cd,Dd] = ssdata(sysd);


%% Form Markov parameters
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


%% Plot comparison of eigenvalues
figure,
subplot(1,3,1)
plot(real(eig(Arf)),imag(eig(Arf)),'bo','markersize',10)
title('Full SVD')
set(gca,'fontsize',16)

subplot(1,3,2)
plot(real(eig(Ar_svds)),imag(eig(Ar_svds)),'rd','markersize',10)
title('SVDS-H')
set(gca,'fontsize',16)

subplot(1,3,3)
plot(real(eig(Ar_rsvdh)),imag(eig(Ar_rsvdh)),'k*','markersize',10)
title('RandSVD-H')
set(gca,'fontsize',16)

