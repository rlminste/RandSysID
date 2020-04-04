function y = blockhankmult(T,x,s,l,m,tflag)
% 
% Computes a matvec with a block Hankel matrix
%
% l = number of outputs 
% m = number of inputs 
% s = number of block rows/columns of block Hankel matrix
% (Hankel matrix is sl x sm)
% 
% tflag = 'transpose' for transpose multiplication, otherwise normal
%
% T is (lm)x(2s-1) matrix representing block Hankel matrix:
%       first s elements is first row r 
%       last s elements is last column c
%       (one overlapped)
%
% x is the vector to be multiplied
%
% y = final matvec product
%
%

% swap l and m for transpose multiplication

if strcmp(tflag,'transpose')  
    temp2 = l;
    l = m;
    m = temp2;
end    

% initialize
y = zeros(s*l,1);
c = zeros(s,1);
r = zeros(1,s);
w = zeros(s,1);
    
    
if strcmp(tflag,'transpose') 
    %transpose multiplication
    for i = 1:l
        for j = 1:m 
            ind = i+(j-1)*l;
            r(1:q) = T(ind,q:-1:1);
            c(1:p) = T(ind,q:end);
            w(1:q) = x(j:m:end);
            Ax = hankelmult(w,c,r);
            y(i:l:end) = y(i:l:end) + Ax;

        end
    end
    
else
    %normal multiplication
    ind = 1;
    for i = 1:l 
        for j = 1:m 
            r(1:q) = T(ind,q:-1:1);
            c(1:p) = T(ind,q:end);
            w(1:q) = x(j:m:end);
            Ax = hankelmult(w,c,r);
            y(i:l:end) = y(i:l:end) + Ax;
            ind = ind +1;
        end
    end
end