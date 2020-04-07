%% Fig8_powersys_markoverr.m
% 
% This code generates Figure 5.8 from the paper
%   'Efficient Randomized Algorithms for
%    Subspace System Identification'
%       -Minster, Saibaba, Kar, Chakrabortty

%load matrices
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

%% Form new markov parameters

markov_rsvdh = cell(1,N);
f = B_rsvdh;
for jj = 1:N
    g = C_rsvdh*f;
    f = A_rsvdh*f;           
    markov_rsvdh{jj} = g;
end

markov_rtera = cell(1,N);
f = B_rtera;
for jj = 1:N
    g = C_rtera*f;
    f = A_rtera*f;           
    markov_rtera{jj} = g;
end

%% Compute relative error in Markov parameters

merr_rsvdh = zeros(1,N);
merr_rtera = zeros(1,N);
for k = 1:N
    merr_rsvdh(k) = norm(markov{k}-markov_rsvdh{k})/norm(markov{k});
    merr_rtera(k) = norm(markov{k}-markov_rtera{k})/norm(markov{k});
end

%% Plot error

figure,
semilogy(1:N,merr_rsvdh,'linewidth',2), hold on
semilogy(1:N,merr_rtera,'k--','linewidth',2)
legend('RandSVD-H','RandTERA')
title('Relative Error in Markov Parameters')
ylabel('Relative Error')
set(gca,'fontsize',18)