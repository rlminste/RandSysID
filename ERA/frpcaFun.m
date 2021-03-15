function [U,S,V] = frpcaFun(Afun,m,n,k,q,mode)
% Afun: function representing product of matrix and transpose
% [m,n] = size(A)
% k = target rank
% q = number of passes through the matrix
% mode = 1 for frPCA, mode = 1 for frPCAt

s = 5;
if mode == 1
    if rem(q,2) == 0
        Q = randn(n,k+s);
        Q = Afun(Q,'notransp');
        if q == 2
            [Q,~,~] = eigSVD(Q);
        else
            [Q,~] = lu(Q);
        end
    else
        Q = randn(m,k+s);
    end
    upper = floor((q-1)/2);
    for i = 1:upper
        if i == upper
            B = Afun(Q,'transp');
            C = Afun(B,'notransp');
            [Q,~,~] = eigSVD(C);
        else
            B = Afun(Q,'transp');
            C = Afun(B,'notransp');
            [Q,~] = lu(C);
        end
    end
    B = Afun(Q,'transp');
    [V,S,U] = eigSVD(B);
    ind = s+1:k+s;
    U = Q*U(:,ind);
    V = V(:,ind);
    S = S(ind);
elseif mode == 2
    if rem(q,2) == 0
        Q = randn(m,k+s);
        Q = Afun(Q,'transp');
        if q == 2
            [Q,~,~] = eigSVD(Q);
        else
            [Q,~] = lu(Q);
        end
    else
        Q = randn(n,k+s);
    end
    upper = floor((q-1)/2);
    for i = 1:upper
        if i == upper
            B = Afun(Q,'notransp');
            C = Afun(B,'transp');
            [Q,~,~] = eigSVD(C);
        else
            B = Afun(Q,'notransp');
            C = Afun(B,'transp');
            [Q,~] = lu(C);
        end
    end
    B = Afun(Q,'notransp');
    [U,S,V] = eigSVD(B);
    ind = s+1:k+s;
    U = U(:,ind);
    V = Q*V(:,ind);
    S = S(ind);
end
end