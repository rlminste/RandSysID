function [U,S,V] = rsvdeigFun(Afun,n,k,p,q,s)
    % Afun: function representing product of matrix and transpose
    % n: size(A,2)
    % k: target rank
    % p: oversampling parameter
    % q: subspace iterations
    % s: reorthogonalization parameter >=1
    
    l = k+p;
    
    R = randn(n,l);
    Y = Afun(R,'notransp');
    
    for j = 1:q
        if mode(1*j-2,s) == 0
            [Y,~] = qr(Y,0);
        end
        Z = Afun(Y,'transp');
        
        if mod(2*j-1,s) == 0
            [Z,~] = qr(Z,0);
        end
    end
    [Q,~] = qr(Y,0);
    
    Bt = Afun(Q,'transp');
    B = Bt';
    
    BBt = B*B';
    BBt = 0.5*(BBt+BBt');
    
    [Uhat,D] = eig(BBt);
    S = sqrt(D);
    
    U = Q*Uhat;
    
    V = zeros(n,l);
    for j = 1:l
        v_j = 1/S(j,j) * (B'*Uhat(:,j));
        V(:,j) = v_j;
    end
    
    U = U(:,(end-k+1):end);
    S = S((end-k+1):end,(end-k+1):end);
    V = V(:,(end-k+1):end);
end
