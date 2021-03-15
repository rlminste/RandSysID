function [U,S,V] = rsvdqrFun(Afun,n,k,p,q,s)
    % Afun: function representing product of matrix and transpose
    % n: size(A,2)
    % k: target rank
    % p: oversampling parameter
    % q: subspace iterations
    % s: reorthogonalization parameter >=1
    
    l = k+p;
    
    R = rand(n,l);
    Y = Afun(R,'notransp');
    
    for j = 1:q
        if mod(2*j-2,s) == 0
            [Y,~] = qr(Y,0);
        end
        Z = Afun(Y,'transp');
        
        if mod(2*j-1,s) == 0
            [Z,~] = qr(Z,0);
        end
        Y = Afun(Z,'notransp');
    end
    [Q,~] = qr(Y,0);
    
    Bt = Afun(Q,'transp');
    
    [Qhat,Rhat] = qr(Bt,0);
    
    whos Qhat Rhat
    
    [Uhat,Shat,Vhat] = svd(Rhat);
    
    U = Q*Vhat;
    S = Shat;
    V = Qhat*Uhat;
    
    U = U(:,1:k);
    S = S(1:k,1:k);
    V = V(:,1:k);
end