function [A,B,C,S] = impulse_era(markov,s,l,m,r,method,parameters)
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
%            - 'rsvdeigH' runs RSVDeig-H
%            - 'rsvdqrH' runs RSVDqr-H
%            - 'frpcaH' runs frPCA-H
%       parameters: structure holding parameters for RandSVD and TERA if using
%            - parameters.rho: oversampling parameters (default 20)
%            - parameters.maxiter: number of subspace iterations (default 1)
%            - parameters.epsilon: threshold for TERA (default 1e-2)
%
% Outputs: 
%       A,B,C: identified system matrices
%       t: timing for the chosen identification method
%
% - Written by Arvind K. Saibaba and Rachel Minster, 2020

if nargin < 7
    %rsvd parameters
    rho = 20;
    maxiter = 1;

    %tera parameters
    epsilon = 1e-2;
else 
    rho = parameters.rho;
    maxiter = parameters.maxiter;
    epsilon = parameters.epsilon;
end
    

switch method
    % Full SVD-ERA
    case 'full'
        H = hankelize(markov,s,l,m);
        size(H)
        [U,S,V] = svd(H,0);
        clear H
    
    % SVDS-H
    case 'svdshankel'
        H = cell2mat(markov);
        T = zeros(l*m,2*s-1);   
        for i = 1:m
            T(i:m:end,:) = H(1:l,i:m:end);
        end
        [U,S,V] = svds(@(x,tflag) svds_blockhank(x,tflag,T,l,m,s),[s*l,s*m],r);
    
    % RandSVD-H    
    case 'randsvdhankel'
        H = markov;
        Hf = makehankelfun(H,p,q,l,m);
        [U,S,V] = rsvdFun(Hf,q*m,r,rho,maxiter);
    
    % TERA
    case 'tera'
        [markov_proj, W1, W2,~] = tera(markov,epsilon);            
        l = size(W1,2); 
        m = size(W2,2);
        H = hankelize(markov_proj,s,l,m);
        [U,S,V] = svd(H,'econ');
    
    % RandTERA    
    case 'randtera'
        N = size(markov,2);
        [markov_proj, W1, W2,~] = tera(markov,epsilon); 
        l = size(W1,2); 
        m = size(W2,2);
        
        H = markov_proj(1:N);
        Hf = makehankelfun(H,s,l,m);
        [U,S,V] = rsvdFun(Hf,s*m,r,rho,maxiter);
    
    % RSVDeig-H
    case 'rsvdeigH'
        H = markov;
        Hf = makehankelfun(H,s,l,m);
        [U,S,V] = rsvdeigFun(Hf,s*m,r,rho,maxiter,1);
    
    % RSVDqr-H
    case 'rsvdqrH'
        H = markov;
        Hf = makehankelfun(H,s,l,m);
        [U,S,V] = rsvdqrFun(Hf,s*m,r,rho,maxiter,1);
    
    % frPCA-H
    case 'frpcaH'
        H = markov;
        Hf = makehankelfun(H,s,l,m);
        [U,S,V] = frpcaFun(Hf,s*l,s*m,r,4,2);
end

 [A,B,C] = era(U,S,V,r,m,l);
 
 if strcmp(method,'randtera')|| strcmp(method,'tera')
     B = B*W2';
     C = W1*C;
 end
end
