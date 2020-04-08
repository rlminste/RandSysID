function [A,B,C,t] = impulse_era(markov,s,l,m,r,method)
%
%runs system identification algorithms
%
% Inputs: 
%       markov: cell array of markov parameters 
%       s: number of block rows/columns of hankel matrix 
%       l: number of outputs
%       m: number of inputs
%       r: target system size (choose size of system n for no model
%       reduction)
%       method: chooses which ID algorithm to run
%            - 'full' runs full SVD-ERA
%            - 'svdshankel' runs SVDS-H
%            - 'randsvdhankel' runs RandSVD-H
%            - 'tera' runs TERA
%            - 'randtera' runs RandTERA
%
% Outputs: 
%       A,B,C: identified system matrices
%       t: timing for the chosen identification method
%
% - Written by Arvind K. Saibaba and Rachel Minster, 2020

%rsvd parameters
rho = 20;
maxiter = 1;

%tera parameters
epsilon = 1e-2;

switch method
    % Full SVD-ERA
    case 'full'
        tic
        H = hankelize(markov,s,l,m);
        size(H)
        [U,S,V] = svd(H,0);
        t = toc;
        clear H
    
    % SVDS-H
    case 'svdshankel'
        tic
        H = cell2mat(markov);
        T = zeros(l*m,2*s-1);   
        for i = 1:m
            T(i:m:end,:) = H(1:l,i:m:end);
        end
        [U,S,V] = svds(@(x,tflag) svds_blockhank(x,tflag,T,l,m,s),[s*l,s*m],r);
        t = toc;
    
    % RandSVD-H    
    case 'randsvdhankel'
        tic
        H = markov;
        Hf = makehankelfun(H,p,q,l,m);
        [U,S,V] = rsvdFun(Hf,q*m,r,rho,maxiter);
        t = toc;
    
    % TERA
    case 'tera'
        [markov_proj, W1, W2,~] = tera(markov,epsilon);            
        l = size(W1,2); 
        m = size(W2,2);
        H = hankelize(markov_proj,s,l,m);
        [U,S,V] = svd(H,'econ');
   
        t=0;
    
    % RandTERA    
    case 'randtera'
        N = size(markov,2);
        [markov_proj, W1, W2,t_newm,~,~] = tera(markov,epsilon); 
        l = size(W1,2); 
        m = size(W2,2);
        
        tic; H = markov_proj(1:N);
        Hf = makehankelfun(H,s,l,m);
        [U,S,V] = rsvdFun(Hf,s*m,r,rho,maxiter);
        t_r = toc;
        t = t_newm+t_r;
end

 [A,B,C] = era(U,S,V,r,m,l);
 
 if strcmp(method,'randtera')|| strcmp(method,'tera')
     B = B*W2';
     C = W1*C;
 end
end
