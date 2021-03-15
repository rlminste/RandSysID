%% Fig9_powersys_response.m
% 
% This code generates Figure 5.9 from the paper
%   'Efficient Algorithms for Eigensystem
%    Realization using Randomized SVD'
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

%% Form identified systems

sys_real = ss(Ad,Bd,Cd,0,ts);
sys_rsvd = ss(Ar_rsvdh,B_rsvdh,C_rsvdh,0,ts);
sys_rtera = ss(Ar_rtera,B_rtera,C_rtera,0,ts);

%% Compute impulse response

t=0:0.007:9;
yt_real = impulse(d2c(sys_real),t);
yt_rsvd = impulse(d2c(sys_rsvd),t);
yt_rtera = impulse(d2c(sys_rtera),t);

%% Plot impulse response

figure
subplot(2,1,1)
plot(t,yt_real(:,2,1),'b','linewidth',2), hold on
plot(t,yt_rsvd(:,2,1),'r--','linewidth',2)
plot(t,yt_rtera(:,2,1),'y:','linewidth',2.5)
grid on
legend('original','RandSVD-H','RandTERA');
ylabel('Rotor speed \omega of Generator 1')
xlabel('time')

subplot(2,1,2)
plot(t,yt_real(:,9,1),'b','linewidth',2), hold on
plot(t,yt_rsvd(:,9,1),'r--','linewidth',2)
plot(t,yt_rtera(:,9,1),'y:','linewidth',2)
grid on
legend('original','RandSVD-H','RandTERA');
ylabel('Rotor speed \omega of Generator 2')
xlabel('time')