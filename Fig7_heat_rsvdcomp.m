%% Fig7_heat_rsvdcomp.m
% 
% This code generates the top left, top right, and bottom right 
% plots in Figure 5.7 from the paper
%   'Efficient Algorithms for Eigensystem
%    Realization using Randomized SVD'
%       -Minster, Saibaba, Kar, Chakrabortty

% load matrices
load steel_A.mat
load steel_B.mat
load steel_C.mat
load steel_E.mat

% Transform system into standard form
L = chol(E,'lower');
A = L\(A/(L'));
B = L\B;
C = C/(L');

% parameters
ts = .007;
m = size(B,2);
l = size(C,1);
r = 20;
s = 1000;
N = 2*s-1;

% convert to discrete
sys = ss(A,B,C,0);
sysd = c2d(sys,ts,'tustin');
[Ad,Bd,Cd,Dd] = ssdata(sysd);

    
%% form markov parameters
markov = cell(1,N);
f = Bd;
for jj = 1:N
    g = Cd*f;
    f = Ad*f;           
    markov{jj} = g;
end

%% System ID algorithms

% SVD-ERA
[Ar,Br,Cr,S] = impulse_era(markov,s,l,m,r,'full');

% RandSVD-H
[Ar_rsvdh,Br_rsvdh,Cr_rsvdh,S_r] = impulse_era(markov,s,l,m,r,'randsvdhankel');

% RSVDeig-H
[Ar_rsvdhv1,Br_rsvdhv1,Cr_rsvdhv1,S_v1] = impulse_era(markov,s,l,m,r,'rsvdeigH');

% RSVDqr-H
[Ar_rsvdhv2,Br_rsvdhv2,Cr_rsvdhv2,S_v2] = impulse_era(markov,s,l,m,r,'rsvdqrH');

% frPCA-H
[Ar_rsvdhfr,Br_rsvdhfr,Cr_rsvdhfr,S_fr] = impulse_era(markov,s,l,m,r,'frpcaH');

sr = diag(S_r);
sv1 = diag(S_v1);
sv2 = diag(S_v2);
sfr = diag(S_fr);

%% Generate new markov parameters

% SVD-ERA
nmarkov = cell(1,N);
A = Ar;
B = Br;
C = Cr;
f = B;
for jj = 1:N
    g = C*f;
    f = A*f;           
    nmarkov{jj} = g;
end

% RandSVD-H
nmarkov_rsvdh = cell(1,N);
A = Ar_rsvdh;
B = Br_rsvdh;
C = Cr_rsvdh;
f = B;
for jj = 1:N
    g = C*f;
    f = A*f;           
    nmarkov_rsvdh{jj} = g;
end

% RSVDeig-H
nmarkov_rsvdhv1 = cell(1,N);
A = Ar_rsvdhv1;
B = Br_rsvdhv1;
C = Cr_rsvdhv1;
f = B;
for jj = 1:N
    g = C*f;
    f = A*f;           
    nmarkov_rsvdhv1{jj} = g;
end

% RSVDqr-H
nmarkov_rsvdhv2 = cell(1,N);
A = Ar_rsvdhv2;
B = Br_rsvdhv2;
C = Cr_rsvdhv2;
f = B;
for jj = 1:N
    g = C*f;
    f = A*f;           
    nmarkov_rsvdhv2{jj} = g;
end

% frPCA-H
nmarkov_rsvdhfr = cell(1,N);
A = Ar_rsvdhfr;
B = Br_rsvdhfr;
C = Cr_rsvdhfr;
f = B;
for jj = 1:N
    g = C*f;
    f = A*f;           
    nmarkov_rsvdhfr{jj} = g;
end

% Compute relative error in the Markov parameters
merr_rsvdh = zeros(1,N);
merr_rsvdhv1 = zeros(1,N);
merr_rsvdhv2 = zeros(1,N);
merr_rsvdhfr = zeros(1,N);
for k = 1:N
    merr_rsvdh(k) = norm(nmarkov{k}-nmarkov_rsvdh{k})/norm(nmarkov{k});
    merr_rsvdhv1(k) = norm(nmarkov{k}-nmarkov_rsvdhv1{k})/norm(nmarkov{k});
    merr_rsvdhv2(k) = norm(nmarkov{k}-nmarkov_rsvdhv2{k})/norm(nmarkov{k});
    merr_rsvdhfr(k) = norm(nmarkov{k}-nmarkov_rsvdhfr{k})/norm(nmarkov{k});
end

%% Plot singular values of H_s
figure,
sfr = sfr(end:-1:1); %frPCA gives an increasing order

subplot(2,2,1)
semilogy(1:20,sr,'linewidth',2), hold on
semilogy(1:20,sv1(end:-1:1),'--','linewidth',2)
semilogy(1:20,sv2,'-.','linewidth',2)
semilogy(1:20,sfr,':','linewidth',2)
xlabel('$j$','interpreter','latex')
ylabel('$\sigma_j$','interpreter','latex')
legend('RandSVD-H','RSVDeig-H','RSVDqr-H','frPCA-H')
title('Singular Values of H_s')
set(gca,'fontsize',18)

%% plot eigenvalues of A

subplot(2,2,2)
plot(real(eig(Ar_rsvdh)),imag(eig(Ar_rsvdh)),'bo','markersize',10), hold on
plot(real(eig(Ar_rsvdhv1)),imag(eig(Ar_rsvdhv1)),'rd','markersize',10)
plot(real(eig(Ar_rsvdhv2)),imag(eig(Ar_rsvdhv2)),'k*','markersize',10)
plot(real(eig(Ar_rsvdhfr)),imag(eig(Ar_rsvdhfr)),'gs','markersize',10)
legend('RandSVD-H','RSVDeig-H','RSVDqr-H','frPCA-H')
xlabel('Re')
ylabel('Im')
xlim([.988 1])
title('Eigenvalues')
set(gca,'fontsize',18)


%% Plot error in Markov parameters
subplot(2,2,4)
semilogy(1:N,merr_rsvdh,'linewidth',2), hold on
semilogy(1:N,merr_rsvdhv1,'--','linewidth',2)
semilogy(1:N,merr_rsvdhv2,'-.','linewidth',2)
semilogy(1:N,merr_rsvdhfr,':','linewidth',2)
title('Markov Parameter Error')
xlabel('$k$','interpreter','latex')
ylabel('Relative Error')
set(gca,'fontsize',18)
