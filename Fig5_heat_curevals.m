%% Fig5_heat_curevals
% 
% This code generates Figure 5.5 from the paper
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

% RandSVD-H
[Ar_rsvdh,Br_rsvdh,Cr_rsvdh,~] = impulse_era(markov,s,l,m,r,'randsvdhankel'); 


%% Plot eigenvalues
figure,
plot(real(eig(Ard1)),imag(eig(Ard1)),'bo','markersize',10), hold on
plot(real(eig(Arf)),imag(eig(Arf)),'rx','markersize',10)
plot(real(eig(Ar_rsvdh)),imag(eig(Ar_rsvdh)),'k+','markersize',10)
xlim([.985 1])
xlabel('Re')
ylabel('Im')
title('CUR-ERA Comparison')
legend('CUR-ERA','SVD','RandSVD-H')
set(gca,'fontsize',18)


