%% Fig6_heat_curmarkoverr
% 
% This code generates Figure 5.6 from the paper
%   'Efficient Algorithms for Eigensystem
%    Realization using Randomized SVD'
%       -Minster, Saibaba, Kar, Chakrabortty
%
% We have used code from the GitHub page associated with the paper
% "System Identifcation via CUR-factored Hankel approximation", 2016, by
% Boris Kramer and Alex A. Gorodetsky 


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


% Parameters
m = size(B,2);     % Inputs
l = size(C,1);     % Outputs
Ns = 500;          % gives s = 1000
Ts = .007;         % Sampling time
r = 20;            % Reduced-order model dimension
s = 1000;
N = 2*s-1;

% convert to discrete
sysfull  = ss(full(A),full(B),full(C),0);
sysfulld = c2d(sysfull,Ts,'tustin');      
[Ad,Bd,Cd,Dd] = ssdata(sysfulld);           


%% generate Markov parameters
markov = cell(1,2*Ns);
markov{1} = Dd;         % First Markov parameter is D term

f = Bd;
for jj = 2:2*Ns         % Computing Markov Parameters iteratively 
    g = Cd*f;
    f = Ad*f;           % markov{jj} = Cd*Ad^(jj-1)*Bd;
    markov{jj} = g;     % + 0.02*norm(markov{ii})*randn(p,m); for noisy data
end

% form block Hankel matrix
H = block_Hankel(markov(2:Ns+1),markov(Ns+1:2*Ns));

%% Low rank matrix approximation through CUR
epsTol = 2e-4;              % Tolerance for maxvol algorithm
rr = r; 
itCUR=1;                    % Since initial column/row selection in maxvol is random,   
tCUR = zeros(itCUR,1);      % can run multiple experiments
Hdiff_full = cell(itCUR,1);
condH = cell(itCUR,1);
Hdiff_iterates = cell(itCUR,1);
HnormIterates = cell(itCUR,1);

relTolDelta = 1e-3; % ||(CUR)_{i} - (CUR)_{i+1} ||_F/||(CUR)_{i}||_F < relTolDelta to stop iteration
compMore =0;
tic
for ii =1:itCUR
    t1 =cputime;
    [I,J,u,s,v, Hdiff_full{ii,jj}, condH{ii,jj}, Hdiff_iterates{ii,jj}, HnormIterates{ii,jj}] = crossApprox(H,rr,relTolDelta,epsTol, compMore);
    t2 = cputime;
    tCUR(ii) = t2-t1;
end


% CUR-based ERA
M1 = H(:,J)*v;
M2 = (u')*H(I,:);
[qO,rO] = qr(M1,0);
[qC,rC] = qr(M2',0);
[u1,s1,v1] = svd(rO*(diag(1./diag(s)))*(rC'));
Uhat = qO*u1;
Vhat = qC*v1;

[Ard1,Brd1,Crd1] = era(Uhat,s1,Vhat,r,m,l); 

%% SVD-ERA and RandSVD-H

% generate markov parameters
markov = cell(1,N);
f = Bd;
for jj = 1:N
    g = Cd*f;
    f = Ad*f;           
    markov{jj} = g;
end


% System ID Algorithms

% Full SVD-ERA
[Arf,Brf,Crf,~] = impulse_era(markov,s,l,m,r,'full');

% SVDS-H
[Ar_svds,Br_svds,Cr_svds,~] = impulse_era(markov,p,no,ni,r,'svdshankel');

% RandSVD-H
[Ar_rsvdh,Br_rsvdh,Cr_rsvdh,~] = impulse_era(markov,s,l,m,r,'randsvdhankel'); 

% generate new markov parameters

% full SVD-ERA
markov_f = cell(1,N);
f = Brf;
for jj = 1:N
    g = Crf*f;
    f = Arf*f;           
    markov_f{jj} = g;
end

% SVDS-H
markov_svds = cell(1,N);
f = Br_svds;
for jj = 1:N
    g = Cr_svds*f;
    f = Ar_svds*f;           
    markov_svds{jj} = g;
end

% RandSVD-H
markov_rsvdh = cell(1,N);
f = Br_rsvdh;
for jj = 1:N
    g = Cr_rsvdh*f;
    f = Ar_rsvdh*f;           
    markov_rsvdh{jj} = g;
end

% CUR-ERA
markov_c = cell(1,N);
f = Brd1;
for jj = 1:N
    g = Crd1*f;
    f = Ard1*f;           
    markov_c{jj} = g;
end

% compute error in markov parameters
err_svds = zeros(1,N);
err_cur = zeros(1,N);
err_rsvdh = zeros(1,N);
for k = 1:N
    err_svds(k) = norm(markov_f{k}-markov_svds{k})/norm(markov_f{k});
    err_rsvdh(k) = norm(markov_f{k}-markov_rsvdh{k})/norm(markov_f{k});
    err_cur(k) = norm(markov_f{k}-markov_c{k})/norm(markov_f{k});
end


%% Plot error in markov parameters
figure,
semilogy(1:N,err_svds,'linewidth',2), hold on
semilogy(1:N,err_rsvdh,'--','linewidth',2)
semilogy(1:N,err_cur,'-.','linewidth',2)
legend('SVDS-H','RandSVD-H','CUR-ERA')
title('Markov Parameter Error')
ylabel('Relative Error')
xlabel('$k$','interpreter','latex')
set(gca,'fontsize',18)

